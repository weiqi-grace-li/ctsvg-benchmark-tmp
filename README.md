# Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework

Code and data for all simulation studies and real-data analyses reported in the paper.

> **Weiqi Li**, **Xinzhou Ge**†, **Yuan Jiang**†
> Department of Statistics, Oregon State University, Corvallis, OR, USA 97331
> † Corresponding authors: [xinzhou.ge@oregonstate.edu](mailto:xinzhou.ge@oregonstate.edu), [yuan.jiang@oregonstate.edu](mailto:yuan.jiang@oregonstate.edu)
>
> Preprint: [bioRxiv 10.1101/2025.11.26.690782](https://www.biorxiv.org/content/10.1101/2025.11.26.690782v1)

## Summary

We benchmark six existing cell-type-specific spatially variable gene (ctSVG) detection methods — C-SIDE, CTSV, spVC, CELINA, STANCE, and MMM — using a realistic, biology-grounded simulation framework built from subcellular-resolution Xenium data, alongside idealized simulations for comparison. A complementary decomposed simulation framework isolates which biological and technical components (cell-type diversity, cell layout, null gene behavior, capture efficiency, spatial pattern realism) drive the sharp drop in method performance from idealized to realistic conditions. These are supplemented by robustness, scalability, and real-data experiments. This repository contains the code for all of these results.

## What this repository contains

Code for every simulation study and real-data analysis in the paper: idealized simulations, realistic (scDesign3-based) simulations, decomposed simulations, robustness experiments, scalability experiments, and real spatial transcriptomics analyses (breast cancer, brain/DLPFC, kidney, lung).

For the realistic simulation, decomposed simulation, and real-data tracks, a **small breast cancer example is included directly in this repository at every pipeline stage** — from preprocessed input, through packaged simulation, to real packaged data — so the full detection pipeline can be run immediately after cloning, without downloading any external data.

## Repository structure

The pipeline for each analysis track runs through up to four stages: **raw data → preprocessing → packaged data → detection**. Where a stage says *(script)*, the code to produce it is included here; datasets beyond the breast cancer example are not committed (large files) and are instead available on Zenodo.

| Track | Raw data → preprocessing → packaged | Detection |
|---|---|---|
| **Idealized** | Simulated directly from assumptions defined in [`simulators/idealized/`](simulators/idealized/) — no external data required | [`run/idealized.R`](run/idealized.R) → `data/detection_results/idealized/` |
| **Realistic** | Download Xenium/Visium data → [`simulators/realistic/preprocessing/`](simulators/realistic/preprocessing/) *(script)* → `data/realistic_preprocessed/` → [`simulators/realistic/realistic_*.R`](simulators/realistic/) *(script)* → `data/realistic_packaged/` | [`run/realistic.R`](run/realistic.R) → `data/detection_results/realistic/` |
| **Decomposed** | Requires a packaged realistic simulation (row above) as input | [`run/decomposed.R`](run/decomposed.R) → `data/detection_results/decomposed/` |
| **Robustness** | Same packaged realistic simulation, with spatial coordinates rotated 0°/30°/60°/90° | Produced as part of [`run/realistic.R`](run/realistic.R); labeled `*_rotate_0`, `*_rotate_30`, `*_rotate_60`, `*_rotate_90` in `data/detection_results/realistic/` |
| **Scalability** | Simulated directly from assumptions defined in [`simulators/scalability/`](simulators/scalability/) — no external data required | [`run/scalability.R`](run/scalability.R) → `data/detection_results/scalability/` |
| **Real data** | Download raw data → [`real_data_processing/*.R`](real_data_processing/) *(script)* → `data/real_packaged/` | [`run/real_data.R`](run/real_data.R) → `data/detection_results/real_data/` |

The breast cancer example is committed at every stage above, so it needs no download. Full-size data for the other tissues/datasets (ovarian, lymph node, lung, DLPFC, kidney) will be available on Zenodo: `[Zenodo DOI — to be added upon publication]`.

### Directory layout

```
.
├── util.R                       shared config, dataset builders, method dispatch, and
│                                 result aggregation used by every run/*.R script
├── method_wrappers/             one R wrapper per benchmarked method
│   ├── CELINA_util.R
│   ├── CSIDE_util.R
│   ├── CTSV_util.R
│   ├── spvc_util.R
│   ├── mmm_util.R
│   ├── stance_util.R
│   └── RCTD_util.R              deconvolution helper used in real-data preprocessing
├── simulators/
│   ├── idealized/                4 idealized simulator scenarios + simulator_setup.csv
│   ├── realistic/                 scDesign3-based simulators, one script per tissue
│   │   ├── preprocessing/          Python: raw Xenium/Visium → preprocessed
│   │   ├── realistic_breast.R
│   │   ├── realistic_lymph.R
│   │   ├── realistic_ovarian.R
│   │   └── realistic_validation.R
│   ├── decomposed/                 scDesign_decompose.R
│   └── scalability/                 scalable_sim.R
├── real_data_processing/         raw real data → packaged, one script per dataset
│   ├── breast.R  lung.R  dlpfc.R  kidney_core.R  kidney_inter.R
├── run/                           orchestration: one script per analysis track
│   ├── idealized.R  realistic.R  realistic_validation.R
│   ├── decomposed.R  scalability.R  real_data.R
├── figures/                       fig1.R … fig6.R — reproduce the paper's figures
├── data/                          all pipeline inputs/outputs (see below)
└── environment/                    requirements.txt (Python)
```

```
data/
├── real_raw/                 raw inputs for real_data_processing/*.R
├── real_packaged/             output of real_data_processing/*.R  (breast example included)
├── realistic_raw/              raw Xenium/Visium bundles for the preprocessing notebooks
├── realistic_preprocessed/     output of the preprocessing notebooks  (breast example included)
├── realistic_packaged/         output of simulators/realistic/realistic_*.R  (breast example included)
└── detection_results/          output of run/*.R
```

## Quick start: run the included example

All four scripts below are **directly runnable immediately after cloning the repo and installing the required packages** (see [Requirements](#requirements)) — no data download needed, since the breast cancer example is already committed.

```bash
git clone https://github.com/weiqi-grace-li/ctsvg-benchmark-tmp.git
cd ctsvg-benchmark-tmp
```

Run any of the following from the repository root (all paths in the code are relative to it):

```bash
Rscript run/idealized.R      # fully synthetic idealized simulations
Rscript run/realistic.R      # breast cancer realistic simulation + robustness rotation sweep
Rscript run/decomposed.R     # decomposed simulation, breast cancer
Rscript run/real_data.R      # real breast cancer spatial transcriptomics data
```

`run/scalability.R` is likewise directly runnable — it is fully synthetic and requires no packaged data.

Results are written to `data/detection_results/<track>/`. Each script is set by default to run only the fastest method, `celina`, for a quick demonstration; uncomment the `run_method` line near the top of the script to run all six benchmarked methods (`celina, stance, spvc, ctsv, cside, mmm`).

## Generating your own realistic simulation from preprocessed data

The breast cancer realistic simulation packaged in `data/realistic_packaged/breast/` (consumed by `run/realistic.R` and `run/decomposed.R` above) was itself produced by [`simulators/realistic/realistic_breast.R`](simulators/realistic/realistic_breast.R), starting from the Python-preprocessed Xenium data already included in `data/realistic_preprocessed/breast/`. Because that preprocessed input is already in the repository, this script can be run directly to reproduce — or regenerate with new random seeds — the packaged simulation:

```bash
Rscript simulators/realistic/realistic_breast.R
```

What it does, step by step:

1. **Assemble** the preprocessed cell-, spot-, and gene-level data into a single working dataset.
2. **Fit `scDesign3`** — a negative-binomial copula model per cell type, over the tissue's spatial coordinates.
3. **Select ground-truth ctSVGs** — genes chosen by cell-type prevalence, spatial autocorrelation (Moran's I), and scDesign3 deviance explained.
4. **Simulate** five new replicate datasets from the fitted model: genes *not* selected as ctSVGs have their spatial signal flattened, so only the selected ctSVGs carry true cell-type-specific spatial patterns; simulated cells are then aggregated into pseudo-spots.
5. **Filter and finalize** — extreme-count outlier genes are removed, and the tissue's spatial triangulation boundary (already included in the repository) is attached, skipping the interactive boundary-picking step the script would otherwise require.

The output — saved to `data/realistic_packaged/breast/` — is exactly the input consumed by `run/realistic.R` and `run/decomposed.R`.

`simulators/realistic/realistic_lymph.R`, `realistic_ovarian.R`, and `realistic_validation.R` follow the same pipeline for other tissues but require their own preprocessed inputs, which are not committed here (see [Extending to other tissues and datasets](#extending-to-other-tissues-and-datasets)).

**Working-directory note:** the R scripts above (`run/*.R`, `simulators/realistic/realistic_*.R`) use paths relative to the **repository root** and must be run from there. The Python preprocessing notebooks in `simulators/realistic/preprocessing/` use paths relative to **their own folder** and must be run from inside that folder instead.

## Extending to other tissues and datasets

Beyond the breast cancer example, the paper's full set of realistic simulations (ovarian, lymph node, Visium-validated breast) and real datasets (lung, DLPFC/brain, kidney core, kidney intermediate) follow the same staged pipeline described in [Repository structure](#repository-structure). Because their raw and intermediate files are large, they are not committed to this repository. To reproduce them yourself:

- **Realistic simulations:** run the relevant notebook in `simulators/realistic/preprocessing/` on the corresponding raw Xenium/Visium data, then the matching `simulators/realistic/realistic_<tissue>.R` script.
- **Real data:** run the corresponding `real_data_processing/<dataset>.R` script on the raw data.

Alternatively, the full-size raw, preprocessed, packaged, and detection-result files for every dataset in the paper will be available on Zenodo: `[Zenodo DOI — to be added upon publication]`.

## Requirements

**R** — no version-pinned environment (e.g. `renv.lock`) is available yet; install the following with your usual method (`install.packages()` / `BiocManager::install()`):

- Core: `Seurat`, `Matrix`, `arrow`, `readr`, `ggplot2`, `dplyr`, `purrr`, `scales`, `Triangulation`, `spatstat`, `patchwork`, `ggforce`
- Simulation: `scDesign3`, `SingleCellExperiment`
- Benchmarked methods: `CELINA`, `STANCE`, `CTSV`, `spVC`, `spacexr` (C-SIDE / RCTD), `MCube` (MMM), `BPST`, `MGLM`, `BiocParallel`
- Real-data preprocessing (Bioconductor): `SpatialExperiment`, `zellkonverter`, `AnnotationDbi`, `org.Hs.eg.db`, `spatialLIBD`, `HDF5Array`

**Python** — only needed for `simulators/realistic/preprocessing/`; see [`environment/requirements.txt`](environment/requirements.txt):
`numpy`, `pandas`, `scipy`, `matplotlib`, `scanpy`, `pyarrow`, `opencv-python`, `ome-types`, `shapely`, `scikit-image`, `tifffile`, `tqdm`.

## License

Released under the [MIT License](LICENSE).
