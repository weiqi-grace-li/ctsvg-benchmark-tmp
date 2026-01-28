# ctsvg-benchmark-tmp (Work in Progress) 
Link to the paper: https://www.biorxiv.org/content/10.1101/2025.11.26.690782v1

Link to the data: (Zenodo - coming soon). The plan is to have `data/` and `result/` folders, where `data/` contains simulated realistic datasets and preprocessed real datasets and `result/` stores ctSVG detection output. 

## Overview  
This repository provides pipeline to reproduce all results in paper **"Benchmarking Cell-Type-Specific Spatially Variable Gene Detection Methods Using a Realistic and Decomposable Simulation Framework"**. We propose an integrated evaluation framework that combines: 
1. Idealized simulations (baseline)
2. Xenium-based realistic simulations (biological complexity)
3. Decomposition-based diagnostic analysis (decompose biological complexity to interpretable components)

## Repository Structure 
This repository is organized into `src/` (core utility functions and class definitions) and numbered execution scripts from 01 to 04. 

> **Note**: The Breast Cancer realistic dataset is pre-packaged in this repository for an immediate **Quick Start**..

<div align="center">

<table>
  <thead>
    <tr>
      <th rowspan="2">Result Type</th>
      <th colspan="2">Data Option: From Scratch</th>
      <th colspan="1">Data Option: From Processed</th>
      <th colspan="1">Result<br></th>
    </tr>
    <tr>
      <th>Raw</th>
      <th>Workflow</th>
      <th>Download</th>
      <th>Workflow</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><b>Idealized</b></td>
      <td>None</td>
      <td>None</td>
      <td>None</td>
      <td><code>03_A</code></td>
    </tr>
    <tr>
      <td><b>Realistic</b></td>
      <td>Xenium</td>
      <td><code>01_A</code>, <code>02_A</code></td>
      <td><code>data/processed/realistic/</code></td>
      <td><code>03_B</code></td>
    </tr>
    <tr>
      <td><b>Decomposed</b></td>
      <td>Xenium</td>
      <td><code>01_A</code>, <code>02_A</code></td>
      <td><code>data/processed/realistic/</code></td>
      <td><code>03_C</code></td>
    </tr>
    <tr>
      <td><b>Real Data</b></td>
      <td>Visium, scRNA</td>
      <td><code>01_B</code></td>
      <td><code>data/processed/real/</code></td>
      <td><code>03_D</code></td>
    </tr>
  </tbody>
</table>

</div>

- **Idealized simulations (`03_A`)**: *no external data required*\
  Generated purely based on user-defined parameters and assumptions.
- **Realistic simulations (`03_B`)**: *require raw source data*\
  Uses real Xenium datasets as raw sources and form realistic datasets through scDesign3 and ground truth ctSVG selection steps. The preprocessing is done via `01_A`, and generation of simulation is done via`02_A`. Users can bypass the preprocessing and generation by directly download simulated realistic datasets from `data/processed/realistic`.  In addition, the repository includes a breast cancer realistic dataset for quick start. 
- **Decomposed simulations (`03_C`)**: *require realistic datasets*\
  Takes realistic datasets as reference and generates decomposed datasets with various combinations of "realness" components (i-vi) based on realistic simulation.
- **Real data (`03_D`)**: *require raw source data*\
  Uses real Visium datasets and reference single-cell datasets used for cell-type deconvolution as raw sources and form preprocessed datasets. The reprocessing is done via `01_B`. Users can bypass the preprocessing by directly download simulated realistic datasets from `data/processed/real/`.

## Usage Scenarios 
This accompanying pipeline is designed with multiple entry points, allowing users to rerun ctSVG detection methods and reproduce results with different levels of data dependency. Users may: 
- (Quick Start) run detection methods on all idealized simulations, and breast cancer realistic simulation and decomposed simulations through breast realistic dataset included in the repository; 
- download processed datasets, run detection methods;
- download raw data from sources, regenerate realistic simulations and preprocessed real data, then run detection methods;
- bypass running detection methods and reproduce figures from precomputed results.
