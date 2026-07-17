# Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework

Code, example data, and analysis workflows accompanying:

> **Preprint:** [bioRxiv 10.1101/2025.11.26.690782](https://www.biorxiv.org/content/10.1101/2025.11.26.690782v1)  
> **Data and code:** [Zenodo archive preview](https://zenodo.org/records/21405877?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6Ijk1OTBiNTRmLTNkYTItNDM2OC05M2U2LWRhN2ZjOWUzMTc1YyIsImRhdGEiOnt9LCJyYW5kb20iOiJmMjljNWRlNzEwNDQ3ZTk1NTEyMjMxYTdhNTRjOTEzYSJ9.qEuq9Z-UuWRAGBkrymTtV3c89YNzCSldUhN-2BukXNCjh26nfGUcbckOibW93ienJZTP3c5uXH9ESZ3rLN-nYQ)

## Overview

This study benchmarks six methods for detecting cell-type-specific spatially variable genes (ct-SVGs):C-SIDE, CTSV, spVC, CELINA, STANCE, MMM.

The evaluation includes a new realistic and decomposable simulation framework, supplemented by robustness, scalability and real-data analyses. The realistic simulations are generated from subcellular-resolution Xenium data using `scDesign3`. The decomposed simulation framework then introduces realistic components individually or in controlled combinations to determine which features—such as cell-type diversity, spatial cell layout, null-gene distributions, capture efficiency, and realistic intra-cell-type spatial patterns—drive performance loss relative to idealized simulations.

**Realistic and decomposed simulations can run with example data included in this repository without downloading additional data (see Quick Start)**. 

## Analysis workflow
The repository supports three levels of reproduction:

1. **Reproduce manuscript figures from archived detection results.**  
   Download detection results from Zenodo and run the scripts in `figures/`.

2. **Rerun ct-SVG detection from packaged datasets.**  
   Download the packaged datasets from Zenodo, or use the included breast-cancer example, and run the scripts in `run/`.

3. **Regenerate packaged datasets from raw data.**  
   Download the original raw data, run the relevant preprocessing scripts, and then generate the packaged simulation or real-data objects.


<p align="center">
  <img src="docs/analysis_pipelines.png" alt="Workflow from raw data through preprocessing, packaged datasets, and ct-SVG detection" width="100%">
</p>



## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/weiqi-grace-li/ctsvg-benchmark-tmp.git
cd ctsvg-benchmark-tmp
```

All analysis scripts assume that they are launched from the repository root.

### 2. Install the required software

The analysis requires R. Python is needed only for preprocessing raw Xenium and Visium files in `simulators/realistic/preprocessing/`.

Restore the exact R environment with:

```bash
Rscript -e 'install.packages("renv"); renv::restore()'
```

Install the Python dependencies with:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r environment/requirements-lock.txt
```

### 3. Run with data included in this repository

#### 3.1 Run ct-SVG detection methods on packaged realistic simulation and real data. 

A small realistic simulation (based on Xenium breast tumor data) and a small processed real-data (Her2+ ST data) are included in the repository. The following commands therefore require no additional data download:

```bash
Rscript run/idealized.R # fully synthetic
Rscript run/realistic.R
Rscript run/decomposed.R
Rscript run/real_data.R
Rscript run/scalability.R # fully synthetic 
```

Outputs are written to:

```text
data/detection_results/<analysis-track>/
```

For a quick smoke test, each run script is configured by default to execute a single comparatively fast method. Edit the `run_method` setting near the beginning of the script to run all six benchmarked methods. Some methods are computationally intensive. See [Computational considerations](#computational-considerations) before launching the complete benchmark.

#### 3.2 Regenerating the realistic breast-cancer simulation

The preprocessed results used for generating breast tumor realistic simulation are also included in the repository. Without downloading any additional data, you can generate the realistic simulation from those preprocessed data. 

```bash
Rscript simulators/realistic/realistic_breast.R
```

using the included preprocessed Xenium files in:

```text
data/realistic_preprocessed/breast/
```

The resulting packaged files are the direct inputs to `run/realistic.R` and `run/decomposed.R`.

## Extending the workflow to the full study

For both realistic simulation and real-data analysis, you can directly download the raw data and run through the preprocessing steps before running ct-SVG detection methods.  Alternatively, the Zenodo archive contain preprocessed inputs, packaged datasets, and detection results. You may begin from any stage shown in the workflow figure rather than repeating all upstream processing.

## Data organization

```text
data/
├── real_raw/
│   Raw inputs for `real_data_processing/*.R`.
│
├── real_packaged/
│   Packaged real-data objects consumed by `run/real_data.R`.
│   The breast-cancer example is included.
│
├── realistic_raw/
│   Raw Xenium/Visium bundles used by the realistic preprocessing scripts.
│
├── realistic_preprocessed/
│   Outputs from `simulators/realistic/preprocessing/`.
│   The breast-cancer example is included.
│
├── realistic_packaged/
│   Outputs from `simulators/realistic/realistic_*.R`.
│   The breast-cancer example is included.
│
└── detection_results/
    Outputs from the scripts in `run/`, or archived results downloaded
    from Zenodo.
```

### Python

Python is required only for the raw-data preprocessing scripts under `simulators/realistic/preprocessing/`.

Dependencies are listed in `environment/requirements.txt` and include:

```text
numpy
pandas
scipy
matplotlib
scanpy
pyarrow
opencv-python
ome-types
shapely
scikit-image
tifffile
tqdm
```


## Citation

When using this repository, please cite:

```bibtex
@article{li2025benchmarking,
  title   = {Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework},
  author  = {Li, Weiqi and Ge, Xinzhou and Jiang, Yuan},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.11.26.690782}
}
```

Please also cite the original ct-SVG methods and source datasets used in the relevant analysis.

## License

This repository is released under the [MIT License](LICENSE).

The original datasets remain subject to the licenses and terms specified by their respective providers.

## Contact

For questions about the code or analyses, please open a GitHub issue.

For correspondence about the study:

- Xinzhou Ge: [xinzhou.ge@oregonstate.edu](mailto:xinzhou.ge@oregonstate.edu)
- Yuan Jiang: [yuan.jiang@oregonstate.edu](mailto:yuan.jiang@oregonstate.edu)
