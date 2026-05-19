# eHMC: Self-Calibrated Randomized Hamiltonian Monte Carlo

This repository contains the full codebase associated with the paper 
**Faster Hamiltonian Monte Carlo by Learning Leapfrog Scale: a self-calibrated randomized solution**
(https://arxiv.org/abs/1810.04449).
It includes:

- **`ehmc`**: main R package that contains the R interface, C++/Rcpp implementation, 
and documentation for the eHMC sampler and its partial-refreshment variant prHMC;
- **`ehmcexamples`**: R/Stan package containing Stan benchmark models;
- **`ehmcproposal`**: Python package for proposal calibration, flow-based approximation, and warmup sample generation;
- experiment scripts for proposal calibration, sampler runs, and result analysis.

---

## Installation guide

### Clone the repository

```bash
git clone https://github.com/jstoehr/eHMC.git
cd eHMC
```

### R packages

From the repository root, install the two R packages. 
The installation requires `devtools`, as well as a working RStan installation 
and C++ toolchain for the `ehmcexamples` package than contains Stan models.

```r
install.packages(c("devtools", "rstan", "rstantools"), dependencies = TRUE)
devtools::install("packages/ehmc", dependencies = TRUE)
devtools::install("packages/ehmcexamples", dependencies = TRUE)
```

If documentation needs to be regenerated for the Stan package, use:

```r
devtools::document(
  "packages/ehmcexamples",
  load_code = rstantools::rstantools_load_code
)
```
For the main package on eHMC, use:
```r
devtools::document("packages/ehmc")
```

### Python package

The proposal calibration code is managed through through the `ehmcproposal`.
The following steps guide you on how to create the conda environment and 
how to install the package and its dependencies.

```bash
cd packages/ehmcproposal
conda env create -f environment.yml
conda activate ehmc
python -m pip install -e .
```

If you need to test the installation:

```bash
python -c "import ehmcproposal; print(ehmcproposal.__file__)"
```

Note that the Python package requires a working R installation because `rpy2` is used to call R and Stan code from Python.
Make sure that R is visible from the same shell.

---

## Experimental workflow

The scripts are written to be launched from the repository root and handle three models
```text
banana
MVNorm
BLP
```

### 1. Calibration of the proposal distribution

Run a single calibration:

```bash
python examples/scripts/run_calibration.py \
  --repo_root . \
  --model banana \
  --seed 1 \
  --n_train 50000 \
  --ess_target 2500 \
  --n_warmup 2000 \
  --resampling 1 \
  --max_tries 100
```

Or use the bash launcher. The launcher accepts command-line overrides. 
Check `examples/scripts/run_calibration.py` for the supported options included.
For example:

```bash
bash examples/scripts/run_calibration.sh \
  --model "banana MVNorm" \
  --seed "1:20" \
  --n_train "50000" \
  --ess_target "2500 5000" \
  --n_warmup "2000 5000" \
  --max_tries 100
```

Calibration creates model-specific output folders under:

```text
examples/<model>/
```
The folder `proposal/` contains the calibrated normalizing flow (proposal distribution) and
checkpoints allowing failed or unfinished calibrations to be resumed:
```text
*_checkpoint.pt
*_flow.pt
```
It also stores summaries of the calibration progression:
```text
*_ess_history.parquet
*_scaled_ess_history.parquet
*_beta_history.parquet
```
Depending on the run, the workflow also creates either `warmup_with_resampling/`
or `warmup_wo_resampling/`, containing files of the form
`*_calib_ehmc_w_*.parquet`. These files store the warmup samples generated
at the end of the PMC procedure and used to construct the empirical distribution
used within eHMC and prHMC.


### 2. Run samplers

The R sampler script supports ehmc, prhmc and nuts. For instance for running, one eHMC sampler manually:

```bash
Rscript examples/scripts/run_sampler.R \
  --repo_root . \
  --model banana \
  --algo ehmc \
  --seed 1 \
  --chains 50 \
  --warmup 2000 \
  --iter 2000 \
  --delta 0.651 \
  --metric diag \
  --n_train 50000 \
  --ess_train 0.8 \
  --ess_target 2500 \
  --resampling 1
```

For running NUTS:

```bash
Rscript examples/scripts/run_sampler.R \
  --repo_root . \
  --model banana \
  --algo nuts \
  --seed 1 \
  --chains 50 \
  --warmup 2000 \
  --iter 2000 \
  --delta 0.651 \
  --metric diag
```

You can also use the bash launcher. The launcher accepts command-line overrides. 
Check `examples/scripts/run_sampler.R` for the supported options included.
For example:

```bash
bash examples/scripts/run_sampler.sh \
  --model "banana MVNorm BLP" \
  --algo "ehmc prhmc nuts" \
  --seed "1:20" \
  --warmup "2000 5000" \
  --delta "0.651 0.8" \
  --n_train 50000 \
  --ess_target "2500 5000"
```

### 3. Analyze experiments

After running the samplers, collect summaries and generate figures with `examples/scripts/analyze_experiments.R`.
The analysis script collects sampler summaries, compares eHMC, prHMC, and NUTS.
It then produces diagnostic plots such as normalized ESJD and normalized effective sample sizes.
Figures are written to `examples/figures/`

---

## Adding a new model

To add a new model, create a folder:

```text
examples/<model>/
```

Each model must provide three components:

```text
examples/<model>/model.stan
examples/<model>/model.R
examples/<model>/model.py
```

### `model.stan`

Defines the target posterior distribution in Stan. 
For a new model, the Stan model must first be compiled and initialized to
create a Stan object used to evaluate the log-density and its gradient:

```r
stan_mod <- rstan::stan_model(
  file = "examples/<model>/model.stan"
)

fit <- rstan::sampling(
  stan_mod,
  data = get_data(),
  chains = 0
)
```
This replace the following code lines in `examples/scripts/run_sampler.R`
```r
fit <- ehmcexamples::get_stanfit(opt$model, data)
```
Then you can define similarly the `target_log_pdf` and `target_grad_log_pdf` functions that will be used by the samplers:
```r
target_log_pdf <- function(pars) {
  rstan::log_prob(fit, pars, FALSE)
}
  
target_grad_log_pdf <- function(pars) {
  rstan::grad_log_prob(fit, pars, FALSE)
}
```

The Python proposal calibration workflow accesses `target_log_pdf` from R 
through `rpy2`. A corresponding
routine should therefore be implemented to expose `target_log_pdf` to Python (in place of
the current `load_r_target_log_pdf`).

### `model.R`

Expected functions:
```r
get_data(data_dir = NULL)
get_pars_name(data = get_data())
```

- `get_data()` returns the model data used to build the Stan target.
- `get_pars_name()` returns parameter names used for outputs and summaries.

### `model.py`

Expected functions:

```python
get_data(data_dir = None, device = None, dtype = None)
get_pars_name(data = None)
rprop_init(n, data = None, seed = None, device = None, dtype = None)
```

- `get_data()` returns the model data in Python format.
- `get_pars_name()` returns parameter names used for exported samples.
- `rprop_init()` samples from the initial proposal distribution and returns
  both samples and their log-density values.

The R and Python implementations should describe the same model and use
consistent parameter ordering.

---

## Citation

If you use this repository, please cite the associated manuscript:

```bibtex
@article{wu:etal:2026,
    title = {Faster Hamiltonian Monte Carlo by Learning Leapfrog Scale}, 
    author = {Changye Wu and Pierre Pudlo and Christian P. Robert and Julien Stoehr},
    year = {2026},
    eprint = {1810.04449},
    archivePrefix = {arXiv},
    primaryClass = {stat.CO},
    url = {https://arxiv.org/abs/1810.04449}, 
    note = {Code available at https://github.com/jstoehr/eHMC}
}
```

---

## License

See the package-specific license files `packages/ehmc/`,
`packages/ehmcexamples/`, `packages/ehmcproposal/`.
