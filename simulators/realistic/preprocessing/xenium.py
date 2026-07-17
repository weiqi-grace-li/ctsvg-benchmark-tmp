import scanpy as sc
from scipy.spatial import cKDTree
from skimage.io import imread as skimread
import pandas as pd
import tifffile
import logging as log
import numpy as np

import importlib
import pseudo_util as pseudo_util
from pathlib import Path

from tqdm import tqdm
from matplotlib.patches import Circle
importlib.reload(pseudo_util)

from pseudo_util import (
    transform_coordinates,
    get_xenium_capture_polygon_um,
    plot_polygons,
    get_visium_capture_polygon,
    generate_space_filling_visium_polygons,
    write_mtx,
)
import matplotlib.pyplot as plt
import gc, os, gzip 
from skimage.measure import points_in_poly
from scipy import sparse 
from scipy.io import mmwrite 
import pyarrow

logger = log.getLogger(__name__)

class Xenium:
    def __init__(self, paths, afline_mat, morph_only = False, selected_genes = None, img = True, cell_type_sheet_name = None):
        '''
        Creates a Xenium object: an aligned Xenium capture ready for pseudo-spot generation.
        :param paths: dictionary with keys:
            1) he_img: path to H&E image; this is the coordinate system and visualization background
            2) cell_info: path to cell_info csv
            3) cell_type: path to a cell-type label file, CSV or Excel (.xlsx/.xls)
            4) cell_count: path to cell x gene count matrix h5
            5) transcript: path to transcript.parquet
            6) morph_img: Xenium morphology image
        :param afline_mat: affine transformation from Xenium coordinates to the H&E image passed above
        :param cell_type_sheet_name: sheet name to read when paths['cell_type'] is an Excel workbook (ignored for CSV)
        :return:
        '''

        self.afline_mat = afline_mat.copy()

        ## immediately load
        self.he_img = skimread(paths['he_img'])
        if not morph_only: 
            self.um_pixel = self.get_resolution_info(paths['he_img'])
        else: 
            self.um_pixel = 1/afline_mat[0,0]
            logger.info(f"directly using morphology conversion {self.um_pixel} um/pixel ...")
        self.transcripts = pd.read_parquet(paths['transcript'])
        self.transcripts["cell_id"] = self.transcripts["cell_id"].apply(
            lambda x: x.decode("utf-8", "ignore") if isinstance(x, (bytes, bytearray)) else x)

        # QC filtering: drop low-quality calls (qv < 20) and transcripts with no assigned cell.
        logger.info(f"Total transcript count: {self.transcripts.shape[0]}...")

        self.transcripts = self.transcripts[(self.transcripts.qv >= 20) & (self.transcripts.cell_id != -1) & (self.transcripts.cell_id != "UNASSIGNED")]

        logger.info(f"After filtering for QC and unassigned cell id, total transcript counts : {self.transcripts.shape[0]}...")

        s = self.transcripts["feature_name"]
        is_bytes = s.map(lambda v: isinstance(v, (bytes, bytearray)))
        self.transcripts.loc[is_bytes, "feature_name"] = s[is_bytes].map(lambda b: b.decode("utf-8"))
        
        if selected_genes is not None:
            selected_genes = {g.decode("utf-8") if isinstance(g, (bytes, bytearray)) else g
                      for g in selected_genes}
            self.transcripts = self.transcripts[self.transcripts["feature_name"].isin(selected_genes)]
            logger.info(f"After filtering for only selected genes, total transcript counts : {self.transcripts.shape[0]}...")
            self.selected_genes = selected_genes

        self.cell_info = pd.read_csv(paths['cell_info']).set_index('cell_id')

        ## add converted cell coordinates
        transcripts_coords_he_image = (
            transform_coordinates(
                self.transcripts[['x_location', 'y_location']].values,
                afline_mat,
            )
        )

        cell_coords_he_image = (
            transform_coordinates(
                self.cell_info[['x_centroid', 'y_centroid']].values,
                afline_mat,
            )
        )

        self.transcripts[["x_transformed", "y_transformed"]] = transcripts_coords_he_image
        self.cell_info[["x_transformed", "y_transformed"]] = cell_coords_he_image

        logger.info("Locations converted.  Now plotting...")
        ## visually examine alignment results
        
        self.xenium_capture_polygon = transform_coordinates(
                get_xenium_capture_polygon_um(paths['morph_img']),
                afline_mat,
            )
        
        if img:
            fig = plt.figure(figsize = (8, 8))
            ax = plt.gca()
            ax.imshow(self.he_img, origin="lower")
            ax.scatter(
                *cell_coords_he_image.T,
                s = 0.1,
                linewidths = 0,
                edgecolors="none",
                marker = "o",
                color = "gold",
            )
            plot_polygons([self.xenium_capture_polygon], ax = ax,
                          facecolor = (0, 0, 0, 0), edgecolor = "black", linewidth = 1)
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.yaxis.set_ticks_position("left")
            ax.xaxis.set_ticks_position("bottom")
            ax.set_aspect("equal", adjustable="box")


        ## save for later, only matters in saving data
        self.cell_type_path = paths['cell_type']
        self.cell_count_path = paths['cell_count']
        self.cell_type_sheet_name = cell_type_sheet_name




    def get_resolution_info(self, path):
        if Path(path).name.lower().endswith(('.ome.tif', '.ome.tiff')):
            with tifffile.TiffFile(path) as tif:
                tags = tif.pages[0].tags
                try:
                    xres = tags["XResolution"].value
                    yres = tags["YResolution"].value
                    unit = tags["ResolutionUnit"].value  # 1=None, 2=inch, 3=cm
                    if unit == 3:
                        x_ratio = 10000/xres[0]*xres[1]
                        y_ratio = 10000 /yres[0] * yres[1]
                    if unit == 2:
                        x_ratio = 2.54 * 10000/xres[0]*xres[1]
                        y_ratio = 2.54 * 10000/yres[0]*yres[1]
                    if x_ratio != y_ratio:
                        logger.warning(f"X ratio  {x_ratio} and Y ratio {y_ratio} mismatch.")

                    logger.info(f"Successfully extracted H&E meta: {x_ratio} um/pixel")

                    return x_ratio

                except KeyError:
                    logger.warning("H&E image doesn't resolution metadata!")
                    return None
        else:
            logger.warning("H&E image doesn't have um to pixel conversion.")
            return None

    def polygons_from_coords(self, coords, diameter, barcodes, img = True):
        ## need to make sure the coords and radius are given in the coordinate of the main h&e
        assert coords.shape[0] == len(barcodes), (f"Expected spot coordinates to match barcodes, got {coords.shape[0]} vs. {len(barcodes)}")

        self.hexagon_polygons = generate_space_filling_visium_polygons(
            spot_coords = coords
        )
        self.generated_capture_polygon = get_visium_capture_polygon(
            spot_centroids = coords,
            spot_diameter= diameter,
        )

        self.spot_names = barcodes
        self.coords = pd.DataFrame(coords, index = self.spot_names, columns = ['x_transformed', 'y_transformed'])


        if img:
            fig = plt.figure(figsize=(8, 8))
            ax = plt.gca()

            ## plotting code
            ax.imshow(self.he_img, origin="lower")
            ax.scatter(
                self.coords['x_transformed'], self.coords['y_transformed'],
                s=1,
                linewidths=0,
                edgecolors="none",
                marker="o",
                color="red",
            )
            plot_polygons(
                self.hexagon_polygons,
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="black",
                linewidth=1,
            )
            plot_polygons(
                [self.generated_capture_polygon],
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="cyan",
                linewidth=1,
            )

            plot_polygons(
                [self.xenium_capture_polygon],
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="black",
                linewidth=1,
            )

            ## plot style adjustments
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.yaxis.set_ticks_position("left")
            ax.xaxis.set_ticks_position("bottom")
            ax.set_aspect("equal", adjustable="box")

    def polygons_from_assumption(self, margin = 5, margin_set = None, diameter = 55, spacing = 100, img = True):
        '''
        Generates a synthetic hexagonal grid of pseudo-spot centers for tissues with no real
        Visium capture to align against. Requires the um/pixel conversion to be known.
        :param margin: inset (in um) from the Xenium capture polygon boundary, used when margin_set is not given
        :param margin_set: explicit [x_min, x_max, y_min, y_max] grid bounds, overrides margin
        :param diameter: spot diameter, in um
        :param spacing: spot-to-spot spacing, in um
        :return:
        '''

        if self.um_pixel is None:
            raise ValueError("um/pixel is None: cannot build assumed grid.")

        spots = []
        dx = spacing/self.um_pixel
        dy = spacing/self.um_pixel * np.sin(np.pi / 3)
        if (margin_set is None): 
            x_min = min(self.xenium_capture_polygon[:, 0]) + margin/self.um_pixel
            x_max = max(self.xenium_capture_polygon[:, 0]) - margin/self.um_pixel
            y_min = min(self.xenium_capture_polygon[:, 1]) + margin/self.um_pixel
            y_max = max(self.xenium_capture_polygon[:, 1]) - margin/self.um_pixel
        else: 
            x_min = margin_set[0]
            x_max = margin_set[1]
            y_min = margin_set[2]
            y_max = margin_set[3]

        row = 0
        y = y_min
        while y <= y_max:
            x_offset = (row % 2) * (dx / 2)
            x = x_min + x_offset
            while x <= x_max:
                if x_min <= x <= x_max and y_min <= y <= y_max:
                    spots.append((x, y))
                x += dx
            y += dy
            row += 1

        self.spot_names = [f"spot{i}" for i in range(1, len(spots) + 1)]

        self.coords = pd.DataFrame(spots, index=self.spot_names, columns=['x_transformed', 'y_transformed'])

        logger.info(f"Generated {self.coords.shape[0]} spots ....")

        self.hexagon_polygons = generate_space_filling_visium_polygons(
            spot_coords = self.coords.values
        )
        self.generated_capture_polygon = get_visium_capture_polygon(
            spot_centroids = self.coords.values,
            spot_diameter= diameter/self.um_pixel,
        )

        if img:
            fig = plt.figure(figsize=(8, 8))
            ax = plt.gca()

            ## plotting code
            ax.imshow(self.he_img, origin="lower")
            ax.scatter(
                self.coords['x_transformed'], self.coords['y_transformed'],
                s=1,
                linewidths=0,
                edgecolors="none",
                marker="o",
                color="red",
            )
            plot_polygons(
                self.hexagon_polygons,
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="black",
                linewidth=1,
            )
            plot_polygons(
                [self.generated_capture_polygon],
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="cyan",
                linewidth=1,
            )

            plot_polygons(
                [self.xenium_capture_polygon],
                ax=ax,
                facecolor=(0, 0, 0, 0),
                edgecolor="black",
                linewidth=1,
            )

            ## plot style adjustments
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.yaxis.set_ticks_position("left")
            ax.xaxis.set_ticks_position("bottom")
            ax.set_aspect("equal", adjustable="box")

    def generate_pseudo_spot(self, radius = 27.5, batch_size = 100000):
        if self.um_pixel is None:
            radius_px = radius
        else:
            radius_px = radius/self.um_pixel

        transcript_array = self.transcripts[["x_transformed", "y_transformed", "feature_name", "transcript_id", "cell_id"]].to_numpy()
        gene_names, gene_indices = np.unique(transcript_array[:, 2], return_inverse=True)


        ## assign transcripts to spots
        tree = cKDTree(self.coords.values)

        num_spot = self.coords.shape[0]
        num_genes = len(gene_names)
        num_transcript = transcript_array.shape[0]

        spot_gene_matrix = np.zeros((num_spot, num_genes), dtype = int)
        assigned_spot = np.empty(num_transcript, dtype = np.int32)
        assigned = np.empty(num_transcript, dtype = bool)

        for i in tqdm(range(0, num_transcript, batch_size), desc = "Processing transcripts", unit = "batch"):
            batch_end = min(i + batch_size, num_transcript)
            batch = transcript_array[i:batch_end]
            distances, spot_indices = tree.query(batch[:, :2])
            valid_mask = distances < radius_px
            valid_spot_indices = spot_indices[valid_mask]
            valid_gene_indices = gene_indices[i:batch_end][valid_mask]

            np.add.at(spot_gene_matrix, (valid_spot_indices, valid_gene_indices), 1)
            assigned_spot[i:batch_end] = spot_indices
            assigned[i:batch_end] = valid_mask

        ## to include information for gene names and spot information, we convert np.array to pd.DF
        self.pseudo_count = pd.DataFrame(spot_gene_matrix, index = self.spot_names, columns = gene_names)
        self.pseudo_count = self.pseudo_count.loc[self.pseudo_count.sum(axis=1) != 0]
        self.pseudo_meta = self.coords.loc[self.pseudo_count.index]
        self.transcripts["assigned_spot"] = self.coords.index[assigned_spot]
        self.transcripts["assigned"] = assigned
        self.radius = radius_px
        logger.info(f"Total transcript {self.transcripts.shape[0]}, "
                 f"in capture transcript {sum(points_in_poly(self.transcripts[['x_transformed', 'y_transformed']].values, self.generated_capture_polygon))},"
                 f"in visium disc transcript {self.pseudo_count.values.sum()}")


    def pseudo_visualize(self, zoom_min = 750, zoom_max = 1000):
        fig, axes = plt.subplots(1, 2, figsize=(16, 8))

        ## first part of the plot

        axes[0].imshow(self.he_img, origin="lower")

        plot_polygons(
            self.hexagon_polygons,
            ax=axes[0],
            facecolor=(0, 0, 0, 0),
            edgecolor="black",
            linewidth=1,
        )
        plot_polygons(
            [self.generated_capture_polygon],
            ax=axes[0],
            facecolor=(0, 0, 0, 0),
            edgecolor="cyan",
            linewidth=1,
        )
        axes[0].scatter(
            self.cell_info["x_transformed"], self.cell_info["y_transformed"],
            s=0.1,
            linewidth=0,
            edgecolor=None,
            marker="o",
            color="gold"
        )
        plot_polygons(
            [self.xenium_capture_polygon],
            ax=axes[0],
            facecolor=(0, 0, 0, 0),
            edgecolor="black",
            linewidth=1,
        )

        axes[0].scatter(self.transcripts.loc[self.transcripts['assigned'], 'x_transformed'],
                        self.transcripts.loc[self.transcripts['assigned'], 'y_transformed'], color="red", s=0.8)

        ## second part of plot
        axes[1].imshow(self.he_img, origin="lower")
        plot_polygons(self.hexagon_polygons, ax=axes[1],
                      facecolor=(0, 0, 0, 0), edgecolor="black", linewidth=1)
        plot_polygons([self.generated_capture_polygon], ax=axes[1],
                      facecolor=(0, 0, 0, 0), edgecolor="cyan", linewidth=1)
        plot_polygons([self.xenium_capture_polygon], ax=axes[1],
                      facecolor=(0, 0, 0, 0), edgecolor="blue", linewidth=1)

        pseudo_spot = self.pseudo_meta[['x_transformed', 'y_transformed']].values
        if (self.radius != None):
            for x, y in pseudo_spot:
                circle = Circle((x, y), radius=self.radius, edgecolor='blue', facecolor='none', linewidth=0.3)
                axes[1].add_patch(circle)
            axes[1].scatter(self.transcripts.loc[self.transcripts['assigned'], 'x_transformed'],
                            self.transcripts.loc[self.transcripts['assigned'], 'y_transformed'], color="red", s=0.8)

        # Set zoom region
        axes[1].set_xlim(zoom_min, zoom_max)  # adjust as needed
        axes[1].set_ylim(zoom_min, zoom_max)  # adjust as needed
        ## plot style adjustments
        for ax in axes:
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.yaxis.set_ticks_position("left")
            ax.xaxis.set_ticks_position("bottom")
            ax.set_aspect("equal", adjustable="box")

        plt.tight_layout()
        plt.show()

    def save_data(self, path, cell_type_id_name = "Barcode", cell_type_cluster_name = "Cluster" ,random_prop = None, sel_spots =None):
        pd.options.mode.copy_on_write = True
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        
        # Select which pseudo spots are "train" (sel_spots) vs "test" (everything else).
        if (random_prop is not None) and (random_prop > 0):
            spots_all = np.asarray(self.pseudo_meta.index, dtype=str)
            k = k = max(1, min(len(spots_all), int(len(spots_all) * float(random_prop))))
            sel_spots = np.sort(np.random.default_rng(42).choice(spots_all, size=k, replace=False))
        elif sel_spots is not None: 
            sel_spots = sel_spots 
        else:
            sel_spots = np.asarray(self.pseudo_meta.index, dtype=str)

        logger.info(f"Selected {sel_spots.shape[0]} out of {self.pseudo_meta.shape[0]} spots...")

        # keep train / test information
        self.pseudo_meta["group"] = "test"
        self.pseudo_meta.loc[sel_spots, "group"] = "train"

        train_cell_ids = self.save_data_transcript_level(path, sel_spots)
        self.save_data_original_cell(path, train_cell_ids, cell_type_id_name, cell_type_cluster_name)
        self.save_data_pseudo_spots(path)

    def save_data_transcript_level(self, path, sel_spots):
        # Keep only the columns needed below.
        tr = self.transcripts.loc[:, ["cell_id", "assigned", "assigned_spot", "feature_name"]].copy()
        for c in ["cell_id", "assigned_spot", "feature_name"]:
            tr[c] = tr[c].astype("category")
        assigned_mask = tr["assigned"].to_numpy()

        # Per-cell-per-spot proportion of assigned transcripts, used to identify train cells
        # (any cell with nonzero mass on a train spot) and exported as cell_spot_prop.
        cell_codes = tr["cell_id"].cat.codes.to_numpy(np.int64)
        spot_codes = tr["assigned_spot"].cat.codes.to_numpy(np.int64)
        n_cells = tr["cell_id"].cat.categories.size
        n_spots = tr["assigned_spot"].cat.categories.size
        
        r = cell_codes[assigned_mask]
        c = spot_codes[assigned_mask]
        cell_spot = sparse.coo_matrix((np.ones(r.size, np.int32), (r, c)), shape=(n_cells, n_spots)).tocsr()
        
        # Per-cell totals (avoid zero)
        cell_total = np.bincount(cell_codes[assigned_mask], minlength=n_cells).astype(np.int32)
        cell_total[cell_total == 0] = 1
        
        # Proportions (float32)
        inv_tot = 1.0 / cell_total.astype(np.float32)
        cell_spot_prop = sparse.diags(inv_tot) @ cell_spot # CSR float32
        
        # Pick training cells (any mass on selected spots)
        spot_idx = pd.Index(tr["assigned_spot"].cat.categories.astype(str))
        sel_spot_idx = spot_idx.get_indexer(sel_spots)        # -1 if unseen
        sel_spot_idx = sel_spot_idx[sel_spot_idx >= 0]
        train_mask = (cell_spot[:, sel_spot_idx].sum(axis=1).A.ravel() > 0)
        train_idx = np.flatnonzero(train_mask)
        
        train_cell_ids = tr["cell_id"].cat.categories[train_idx].astype(str)
        logger.info(f"selected total {len(train_cell_ids)} cells...")
        sel_spot_ids = tr["assigned_spot"].cat.categories[sel_spot_idx].astype(str)
        
        # Write Matrix Market + id lists
        
        write_mtx(path/"cell_spot_prop",
                  cell_spot_prop[train_idx][:, sel_spot_idx],
                  train_cell_ids, sel_spot_ids)
        
        logger.info("Finished writing cell spot proportion...")
        del cell_spot, cell_spot_prop; gc.collect()

        # Per-cell-per-gene dropout fraction: share of each gene's transcripts in a cell
        # that were not assigned to any pseudo spot.
        gene_codes = tr["feature_name"].cat.codes.to_numpy(np.int64)
        n_genes = tr["feature_name"].cat.categories.size

        M_total = sparse.coo_matrix((np.ones(cell_codes.size, np.int32), (cell_codes, gene_codes)),
                                    shape=(n_cells, n_genes)).tocsr().astype(np.float32)

        M_drop  = sparse.coo_matrix((np.ones((~assigned_mask).sum(), np.int32),
                                     (cell_codes[~assigned_mask], gene_codes[~assigned_mask])),
                                    shape=(n_cells, n_genes)).tocsr().astype(np.float32)
        M_drop_frac = M_drop.multiply(1.0)
        M_total_inv = M_total.copy()
        M_total_inv.data = 1.0 / M_total_inv.data
        M_drop_frac = M_drop_frac.multiply(M_total_inv)

        genes = tr["feature_name"].cat.categories.astype(str)
        write_mtx(path/"gene_dropout_frac", M_drop_frac[train_idx], train_cell_ids, genes)
        logger.info("Finished writing cell-gene dropout fraction...")

        del M_drop, M_drop_frac

        return train_cell_ids

    def save_data_original_cell(self, path, train_cell_ids, cell_type_id_name, cell_type_cluster_name):
        if Path(self.cell_type_path).suffix.lower() in (".xlsx", ".xls"):
            cell_type = pd.read_excel(
                self.cell_type_path,
                sheet_name=self.cell_type_sheet_name,
                usecols=[cell_type_id_name, cell_type_cluster_name],
            )
        else:
            cell_type = pd.read_csv(self.cell_type_path, usecols=[cell_type_id_name, cell_type_cluster_name])
        cell_type = (cell_type.set_index(cell_type_id_name, verify_integrity=True))
        cell_type.index = cell_type.index.astype(str)

        cell_info = self.cell_info[["x_transformed", "y_transformed"]].copy()
        cell_info.index = cell_info.index.astype(str)

        cell_meta = cell_info.join(cell_type, how="left")
        cell_meta = cell_meta.loc[train_cell_ids, ["x_transformed","y_transformed",cell_type_cluster_name]]
        cell_meta.columns = ["x", "y", "type"]
        cell_meta.to_parquet(path/"cell_meta.parquet", compression="zstd")
        del cell_meta, cell_type; gc.collect()
        
        ad = sc.read_10x_h5(self.cell_count_path)
        X = ad.X.tocsr()
        mask = np.isin(ad.obs_names.astype(str), train_cell_ids)

        outd = path / "cell_counts_10x"
        outd.mkdir(exist_ok=True)
        
        with gzip.open(outd / "matrix.mtx.gz", "wb") as f:
            mmwrite(f, X[mask].T.tocoo())  # genes x cells
        pd.Series(ad.obs_names[mask].astype(str)).to_csv(
            outd / "barcodes.tsv.gz", index=False, header=False, compression="gzip"
        )
        
        pd.DataFrame({
            0: ad.var_names.astype(str),           # feature_id
            1: ad.var_names.astype(str),           # feature_name
            2: "Gene Expression"                   # feature_type
        }).to_csv(outd / "features.tsv.gz", sep="\t", index=False, header=False, compression="gzip")
        
        del ad, X; gc.collect()
        
        logger.info("Finished saving training cell meta and count information...")
        
        
    def save_data_pseudo_spots(self, path):
        self.pseudo_count.to_parquet(path/"pseudo_count.parquet", compression="zstd")
        self.pseudo_meta.to_parquet(path/"pseudo_meta.parquet", compression="zstd")
        self.pseudo_meta.to_csv(path/"pseudo_meta.csv")
        logger.info("Finished saving pseudo meta and count information")
        
        
        
        
        
        
        
        
        
        

        

