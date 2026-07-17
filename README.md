# ctSVG benchmark — repository structure

This repo benchmarks six spatially-variable-gene (SVG) detection methods — CELINA, spVC,
C-SIDE, MMM, CTSV, STANCE — across idealized simulations, scDesign3-based realistic
simulations (Xenium-derived breast/ovarian/lymph, plus a Visium-validated breast dataset),
decomposed simulations, scalability tests, and five real datasets (breast, DLPFC/brain, kidney
core, kidney inter, lung).

This repo is self-contained: every `source()` / `read.csv()` / `load()` / `save()` path inside
it resolves to somewhere inside this folder, assuming scripts are run with this folder itself
as the working directory. There are no references to sibling folders outside this repo.

## Layout

```
util.R                per-simulation-type config, dataset builders, method dispatch
                        (run_all_tests()), and evaluation/aggregation (collect_details(),
                        analyze_result(), write_sheet()) — one file, deliberately not split up

method_wrappers/      per-method R wrappers called by util.R's run_all_tests():
                        CELINA_util.R  CSIDE_util.R  CTSV_util.R  spvc_util.R
                        mmm_util.R  stance_util.R  (the 6 benchmarked methods)
                        RCTD_util.R  (deconvolution helper used in real-data prep)

simulators/
  idealized/           the 4 idealized simulator scenario classes + simulator_setup.csv
  realistic/           realistic_breast.R  realistic_ovarian.R  realistic_lymph.R
                        realistic_validation.R  (scDesign3 fitting + ctSVG gene selection,
                        one script per tissue)  pseudo_spot_simulator.R
    preprocessing/      Python: xenium.py, companion_functions.py, pseudo_util.py,
                        preprocess_realistic.ipynb, preprocess_realistic_validation.ipynb
  decomposed/          scDesign_decompose.R
  scalability/         scalable_sim.R

real_data/             breast.R  dlpfc.R  kidney_core.R  kidney_inter.R  lung.R
                        (dlpfc.R produces the "brain" results; kidney_core/kidney_inter
                        share one raw input but package + report separately)

run/                   execution/orchestration scripts, one per simulation category:
                        idealized.R  realistic.R  realistic_validation.R  decomposed.R
                        scalability.R  real_data.R

figures/               fig1.R ... fig6.R, reading from data/ and writing to figures/output/
                        + the small aggregated CSVs fig4.R reads directly: decomposition_dim.csv

data/                  see "The data/ folder" below

environment/           requirements.txt (Python, package names only)
```

## The `data/` folder

`data/` holds everything the pipeline reads and writes, organized by stage:

```
data/
  real_raw/                  raw inputs for real_data/*.R (not committed — large; populate
                              yourself before running real_data/*.R)
    breast/  kidney/  lung/  (kidney/ is shared by kidney_core.R and kidney_inter.R)

  real_packaged/              output of real_data/*.R: filter1.RData, filter2.RData (where
                              applicable), Tr.cell.RData, per dataset
    lung/  kidney_core/  kidney_inter/  breast/  dlpfc/

  realistic_raw/              raw Xenium/Visium bundles for the preprocessing notebooks
                              (not committed — large; populate yourself before running
                              simulators/realistic/preprocessing/*.ipynb)
    breast/  visium/  xenium-rep1/  xenium-rep2/  Cell_Barcode_Type_Matrices.xlsx
    ovarian/                   single Xenium bundle (no Visium/replicate-2/cell-type excel)
    lymph/                     single Xenium bundle
    fig3/                       breast_cells.csv, ovarian_cells.csv, lymph_cells.csv — the
                              per-cell coordinate + cell-type tables fig3.R's cell-type
                              scatter panels read directly

  realistic_preprocessed/     output of simulators/realistic/preprocessing/*.ipynb, input to
                              simulators/realistic/realistic_*.R
    breast/  ovarian/  lymph/  validation/

  realistic_packaged/         output of simulators/realistic/realistic_*.R: breast_small_ori /
                              lymph_small_ori / etc., breast_boundary.RData, data_new_v3_<seed>.RData
                              (5 seeds per tissue), extreme_genes_over300.RData, special_genes.RData
    breast/  ovarian/  lymph/  validation/

  detection_results/          output of run/*.R (what fig2.R/fig3.R/fig5.R read back)
    idealized/  realistic/  decomposed/  scalability/  realistic_validation/  real_data/
```

`real_raw/`, `realistic_raw/{breast,ovarian,lymph}` (excluding `fig3/`), and most of
`realistic_packaged/`/`realistic_preprocessed/` beyond `breast/` are scaffolded as empty
directories — they hold large data that isn't committed here. Run the corresponding
`simulators/realistic/preprocessing/*.ipynb` notebooks and `simulators/realistic/realistic_*.R`
/ `real_data/*.R` scripts to populate them; `detection_results/` fills in as you run `run/*.R`.

## Not yet added

- **`LICENSE`** — no license file has been created; choose one and add it before publishing.
- **R dependency lockfile** — `environment/requirements.txt` covers the Python preprocessing
  imports (package names only, unpinned). No `renv.lock` has been generated, since that
  requires actually running the R pipeline in a controlled environment to capture real package
  versions. Recommend running `renv::init()` / `renv::snapshot()` once the pipeline has been
  run end-to-end at least once.
