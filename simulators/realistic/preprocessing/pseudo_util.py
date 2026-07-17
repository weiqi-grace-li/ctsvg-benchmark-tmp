import math
from collections import namedtuple, Counter

import numpy as np
import scanpy as sc
import pandas as pd

## Image manipulation and geometry
from cv2 import perspectiveTransform
from ome_types import from_tiff
from skimage.measure import points_in_poly

from scipy.spatial import KDTree, ConvexHull
from scipy.sparse import csc_matrix
from shapely.geometry import LinearRing

## Plotting imports
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, to_hex, Colormap
from matplotlib.axes import Axes
import gzip 
from scipy import sparse 
from scipy.io import mmwrite 


def transform_coordinates(
    coords: np.ndarray[np.float32, (None, 2)], transform_matrix: np.ndarray[np.float32]
) -> np.ndarray[np.float32, (None, 2)]:
    """Transforms coordinates using transform.
    https://en.wikipedia.org/wiki/Homography_(computer_vision)
    """
    return perspectiveTransform(
        coords.reshape(-1, 1, 2).astype(np.float32), transform_matrix
    )[:, 0, :]


def get_xenium_capture_polygon_um(
    path_to_xenium_morph: str,
) -> np.ndarray[int, (4, 2)]:
    """Produces a polygon which describes the Xenium morphology image shape.
    Polygons are easily transformable. Thus, this polygon is useful for sub-setting coordinates
    that lie within an image under arbitrary transformation, and of course for debug purposes.
    """
    morphology_metadata = from_tiff(path_to_xenium_morph)
    image_shape_px = np.array(
        [
            morphology_metadata.images[0].pixels.size_x,
            morphology_metadata.images[0].pixels.size_y,
        ]
    )
    image_shape_um = (
        image_shape_px * morphology_metadata.images[0].pixels.physical_size_x
    )
    xenium_capture_polygon_um = np.array(
        [(0, 0), (0, image_shape_um[1]), image_shape_um, (image_shape_um[0], 0)]
    )
    return xenium_capture_polygon_um


def plot_polygons(
    poly_coords: list[list[float]],
    ax: Axes = None,
    auto_set_lim: bool = True,
    **patch_collection_kwargs
):
    """
    Plots a set of polygons. Polygons may be styled using patch_collection_kwargs passed to Collection.
    e.g.
    array, facecolor, edgecolor, linewidth....
    https://matplotlib.org/stable/api/collections_api.html#matplotlib.collections.Collection

    pitfall:
    when auto_set_lim == False, a plot will no rescale to include polygons by default.
    """
    if ax is None:
        ax = plt.gca()

    ax.add_collection(
        polygons_coords_to_patch_collection(poly_coords, **patch_collection_kwargs)
    )
    if auto_set_lim:
        ax.autoscale_view()

def polygons_coords_to_patch_collection(
    poly_coords: list[list[float]], **patch_collection_kwargs
) -> PatchCollection:
    """
    Generates a matplotlib PatchCollection from a set of polygons. A helper function for plot_polygons.
    """
    patches = [Polygon(hexagon) for hexagon in poly_coords]
    collection = PatchCollection(patches, **patch_collection_kwargs)
    return collection


def get_median_spot_to_spot_distance_from_centroids(
    centroids: np.ndarray[np.float32, (None, 2)]
) -> float:
    """
    Calculates the median distance between all coordinates and their closest neighbor.
    Usefull for quickly determining the median spot <-> spot distance for Visium spots in the current coordinate system.
    Can also be directly calculated using the full resolution image pixel size and the theoretical spot <-> spot distance.
    This distance is used to generate a rough estimate of the Visium capture area and to produce space filling
    visium spot Hexagons.
    """
    tree = KDTree(centroids)
    distances, neighbor = tree.query(centroids, k=2)
    return np.median(distances[0, 1])

def hex_corner_offset(corner: int, size: float) -> np.ndarray[float, (2,)]:
    """
    Helper function for polygon_corners, calculates the offset of each hexagon vertex from the center.
    """
    angle_deg = 60 * corner - 30
    angle_rad = math.pi / 180 * angle_deg
    offset = np.array([size[0] * math.cos(angle_rad), size[1] * math.sin(angle_rad)])
    return offset


def polygon_corners(
    center: np.ndarray[float, (2,)], size: float
) -> np.ndarray[float, (6, 2)]:
    """
    Generates the polygon vertices for a hexegon of size at origin center.
    Used for plotting Visium space filling Hexagons.
    Reference for working with Hexagons: https://www.redblobgames.com/grids/hexagons/
    Visium spots are organized such that space filling hexagons are POINTY.
    """
    return np.array([center + hex_corner_offset(i, size) for i in range(0, 6)])


def generate_space_filling_visium_polygons(
    spot_coords: np.ndarray[np.float32, (None, 2)], spot_to_spot_distance: float = None
) -> list[np.ndarray[float, (6, 2)]]:
    """
    Generates hexagon polygons for every visium spot, by default the polygons are fully space filling.
    That is to say the size of the hexagon is set using the spot <-> spot distance.
    These polygons are for display purposes only and should not be used for direct analysis.
    """
    if spot_to_spot_distance is None:
        spot_to_spot_distance = get_median_spot_to_spot_distance_from_centroids(
            spot_coords
        )
    size = spot_to_spot_distance / math.sqrt(3)
    visium_hexegon_polygons = np.array(
        [polygon_corners(coord, [size, size]) for coord in spot_coords]
    )
    return visium_hexegon_polygons

__OUTSIDE_VISIUM_CAPTURE_AREA_BARCODE__ = "Outside_Visium_Capture_Area"

def bin_xenium_data_to_visium_spots(
    visium_spot_centroids: np.ndarray[np.float32, (None, 2)],
    visium_spot_diameter: float,
    visium_barcodes: np.ndarray[str, (None, 2)],
    xenium_coords: np.ndarray[np.float32, (None, 2)],
    xenium_to_visium_transform: np.ndarray[np.float32],
) -> np.ndarray[str]:
    """
    Assigns xenium coordinates to visium spots using a transform into the visium coordinate system.
    Coordinates are assigned by nearest neighbor.
    e.g. each Xenium transcript/cell etc is assigned to the closest Visium barcode after transformation.
    Also filters coordinates that lie outside an estimated capture area.
    """
    xenium_coords_in_visium = transform_coordinates(
        xenium_coords, xenium_to_visium_transform
    )
    kdtree = KDTree(visium_spot_centroids)
    _, neighbors = kdtree.query(xenium_coords_in_visium, k=1)
    barcode_assignments = visium_barcodes[neighbors].astype(object)

    visium_capture_polygon_full_res = get_visium_capture_polygon(
        spot_centroids=visium_spot_centroids, spot_diameter=visium_spot_diameter
    )
    barcode_assignments[
        ~points_in_poly(xenium_coords_in_visium, visium_capture_polygon_full_res)
    ] = __OUTSIDE_VISIUM_CAPTURE_AREA_BARCODE__
    return np.array(barcode_assignments, dtype = object)

def get_visium_capture_polygon(
    spot_centroids: np.ndarray[np.float32, (None, 2)], spot_diameter: float
) -> np.ndarray[np.float32, (None, 2)]:
    """
    Generates a polygon representing the visium capture area; It calculates the convex hull of all Visium spots
    then pads the convex hall by the diameter of one spot. This is used for filtering data outside the capture area
    as well as visualization.

    https://en.wikipedia.org/wiki/Convex_hull
    """
    convex_hull = ConvexHull(spot_centroids)
    convex_poly = convex_hull.points[convex_hull.vertices]

    padded_convex_polly = np.array(
        LinearRing(convex_poly).buffer(spot_diameter).exterior.coords.xy
    ).T
    return padded_convex_polly


def unique_encode(
    array: np.ndarray[str], return_inverse: bool = True
) -> (np.ndarray[str], np.ndarray[int]):
    """np.unique is slower then list(set(array)) thus this is an emulator function of np.unique.
    for fast encoding of categorical arrays.
    """
    unique_values = np.array(list(set(array)))
    if not return_inverse:
        return unique_values
    inverse_map = {val: i for i, val in enumerate(unique_values)}
    inverse_codes = np.array([inverse_map[i] for i in array])
    return (unique_values, inverse_codes)

def generate_dataframe_from_transcript_assignments(
    barcodes: np.ndarray[str],
    feature_names: np.ndarray[str],
) -> sc.AnnData:
    """Generates a count matrix as an anndata object from a set of barcodes and feature names.
    e.g. given transcript cell assignments and gene names it will produce a barcode feature matrix.
    """
    valid_feature_names, feature_codes = unique_encode(
        feature_names, return_inverse=True
    )
    valid_barcodes, barcode_codes = unique_encode(barcodes, return_inverse=True)

    mat_dict = Counter(zip(barcode_codes, feature_codes))
    df = pd.DataFrame(
        0, index=valid_barcodes, columns=valid_feature_names, dtype=np.int32
    )
    for (bc, feat), count in mat_dict.items():
        df.iat[bc, feat] = count

    return df

def write_mtx(base, S, row_ids, col_ids):
    S = S.tocoo().astype(np.float32)
    with gzip.open(f"{base}.mtx.gz", "wb") as f: 
        mmwrite(f, S, field="real")
    pd.Series(row_ids).to_csv(f"{base}.rows.tsv.gz", index=False, header=False, compression="gzip")
    pd.Series(col_ids).to_csv(f"{base}.cols.tsv.gz", index=False, header=False, compression="gzip")    
