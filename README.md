# Sequential Anomaly Detection with Online Missed Discovery Rate Control

Reference implementation and reproduction code for the paper

> **Sequential Anomaly Detection with Online Missed Discovery Rate Control**
> Yang Tian, Wenguang Sun, Bowen Gang

The repository implements **OMDRC** — an online decision rule that screens a
data stream one observation at a time while keeping the *missed discovery rate*
(MDR) below a user-specified level $\alpha$ at every time point. Control is
achieved through an "earn-and-spend" capacity (alpha-wealth) recursion: each
rejection earns capacity $\alpha\,\mathrm{Lmdr}_t$, each non-rejection spends
$(1-\alpha)\,\mathrm{Lmdr}_t$, and a non-rejection is only allowed when the
remaining capacity $C_t$ can absorb its cost **and** the current score passes
the causal local no-free-riding barrier in Algorithms 1 and 2.

---

## What is in this repository

| | |
|---|---|
| `code-semi/OMDRC.R`, `application/OMDRC.R` | The core library (identical copies). Every experiment sources one of them. |
| `code-semi/` | All simulation experiments (main text + appendix). |
| `application/CCFD/` | The real-data experiment (credit-card fraud detection). |
| `appendix-fdr-comparison/` | Appendix comparisons against FDR-based procedures. |
| `figures/` | The final figures as they appear in the manuscript (PDF). |
| `env/` | Scripts to install and audit the exact R environment used for the paper. |

**What is deliberately *not* in this repository**

* **Raw data.** The credit-card dataset is redistributed under a licence that
  does not permit mirroring, and the file (144 MB) exceeds GitHub's 100 MB
  per-file limit. See [Data](#data).
* **`*.rds` result caches** (~900 MB in the authors' working tree). Every script
  regenerates its own cache on first run.
* Exploratory probes, diagnostic one-offs, run logs, and real-data experiments
  that did not pass their preregistered validation gates. Only the code behind
  the published figures is included here.

---

## Installation

The paper's results were produced with **R 4.1.2** and the package versions below.

```
mgcv       1.8.39     data.table  1.14.2     foreach       1.5.2
glmnet     4.1.3      ggplot2     3.3.5      doParallel    1.0.17
onlineFDR  2.2.0      dplyr       1.0.8      iterators     1.0.14
REBayes    2.60       tidyr       1.2.0      RColorBrewer  1.1.2
isotree    0.5.22     patchwork   1.1.1      scales        1.1.1
kedd       1.0.3      ggpubr      0.4.0      Matrix        (base R)
```

On a Debian/Ubuntu machine, `env/setup_r_environment.sh` installs all of the
above (it must run as root):

```bash
sudo bash env/setup_r_environment.sh     # apt packages + pinned source builds
Rscript env/audit_r_environment.R        # verify: parses every script, smoke-tests each package
```

`REBayes 2.60`, `isotree 0.5.22` and `onlineFDR 2.2.0` are installed from source
tarballs, which the setup script expects in `env/sources/` and verifies against
`env/checksums.md5`. The tarballs are **not** shipped here; download them from
the CRAN archive:

```
https://cran.r-project.org/src/contrib/Archive/REBayes/REBayes_2.60.tar.gz
https://cran.r-project.org/src/contrib/Archive/isotree/isotree_0.5.22.tar.gz
https://cran.r-project.org/src/contrib/Archive/onlineFDR/onlineFDR_2.2.0.tar.gz
```

If you only want to run the code (not reproduce the exact environment),
`install.packages()` on recent versions of these packages is sufficient.

---

## Quick start

Every script resolves paths from a single environment variable, `OMDRC_ROOT`,
which must point at the root of this repository. If it is unset, the current
working directory is used instead.

```bash
export OMDRC_ROOT="/path/to/this/repository"

# Main-text simulation, Setting (a)
Rscript "code-semi/nofreeride-setting1,2,3/nfr_setting1_fig.R"
```

A minimal use of the library on your own data:

```r
source(file.path(Sys.getenv("OMDRC_ROOT"), "code-semi", "OMDRC.R"))

# z0, z1   labeled null / non-null reference samples
# z_ini    initial unlabeled warm-up batch (size K0)
# z        the incoming stream
res <- OMDRC_DD(z = z, z_ini = z_ini, z0 = z0, z1 = z1,
                alpha = 0.10, ratio_method = "gam",
                pi_bounds = c(0.01, 0.99),
                D_mode = "growing", D_beta = 0.6, D_min = 10,
                w = 100)

res$de        # 1 = rejected (flagged as a discovery), 0 = not rejected
res$Lmdr      # estimated local MDR of each streamed observation
res$pi_hat    # online prevalence estimate pi_t
res$capacity  # capacity (alpha-wealth) path C_t

# Oracle benchmark, when the true local MDR values are known
OMDRC_OR(x.Lmdr = true_Lmdr, x.Lmdr_ini = true_Lmdr_ini,
         alpha = 0.10, w = 100)$de
```

`ratio_method` selects the density-ratio estimator (`"gam"`, `"kde"`,
`"ulsif"`, `"classifier"`); extra arguments for it go through `ratio_control`.
The growing window is controlled by `D_mode`, `D_beta` and `D_min`, and the
prevalence estimate is clipped to `pi_bounds`.

---

## Reproducing the figures

Scripts write their output next to themselves. `figures/` holds the versions
that were compiled into the manuscript, so you can diff your run against them.

| Manuscript figure (`\label`) | Graphic file | Script(s) to run, in order |
|---|---|---|
| `fig:sim-semi` — Settings (a)–(c), MDR/FDR over time | `NFR_Setting123_3x2.pdf` | `code-semi/nofreeride-setting1,2,3/nfr_setting1_fig.R`, `nfr_setting2_fig.R`, `nfr_setting3_fig.R`, then `combine_nfr_settings_3x2.R` |
| `fig:capacity` — dynamic budget allocation mechanism | `Figure2_oracle_mechanism_setting2_surge_pi_NFR.pdf` | `code-semi/Figure2/Figure2_nfr.R` |
| `fig:sample_size_impact` — sensitivity to $n$, $\beta$, $K_0$, $w$ | `Rplot_NFR_D_n_sensitivity.pdf` | `code-semi/sensity/nfr_D_n_vary.R` |
| `fig:fraud_detection` — credit-card fraud application | `Figure3_CCFD_online_comparison_revised_singlecolumn.pdf` | `application/CCFD/CCFD.R` (requires the dataset) |
| density-ratio estimator ablation (KDE / uLSIF / PC-DRE) | `NFR_samplesize_kde_ulsif_gam_2x2_timeline.pdf` | `code-semi/nofreeride-setting1,2,3/nfr_samplesize_kde_ulsif_gam.R` |
| `fig:highdim_ratio` — high-dimensional DRE strategies | `compare_gam_dim_timeline.pdf` | `code-semi/nofreeride-setting1,2,3/nfr_compare_gam_dim.R` |
| `fig:setting1_varied_k` | `NFR_Setting1_varied_k.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting1_varied_k.R` |
| `fig:setting1_varied_pi` | `NFR_Setting1_varied_pi.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting1_varied_pi.R` |
| `fig:setting2_varied_mu` | `NFR_Setting2_varied_mu.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting2_varied_mu.R` |
| `fig:setting2_varied_pi` | `NFR_Setting2_varied_pi.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting2_varied_pi.R` |
| `fig:setting3_varied_d` | `NFR_Setting3_varied_d.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting3_varied_d.R` |
| `fig:setting3_varied_mu` | `NFR_Setting3_varied_mu.pdf` | `code-semi/setting1,2,3-appendix/nfr_setting3_varied_mu.R` |
| `fig:fdr_mdr_tradeoff` — why FDR optimality ≠ MDR control | `FDR_MDR_trade_off.pdf` | `appendix-fdr-comparison/MDR_FDR_trade_off.R` |
| `fig:offline_adadetect_mdr_comparison` — MDRC vs AdaDetect | `OfflineMDRC_vs_AdaDetect_2x2_FDR_left_MDR_right.pdf` | `appendix-fdr-comparison/R1_AdaDetect_OfflineMDRC.R` |

### The no-free-riding barrier (Algorithms 1 and 2)

A capacity-only rule can let an observation with a large `Lmdr` be missed simply
because enough capacity had accumulated earlier. The local barrier removes this
"free-riding" effect while preserving $C_t \ge 0$, so MDR validity is unaffected.

The public `OMDRC_OR()` and `OMDRC_DD()` functions apply this barrier directly.
The manuscript reproduction scripts use explicit `nfr_` names and write
`NFR_`-prefixed outputs so they remain distinguishable from archived pre-revision
experiments:

```
code-semi/nofreeride-setting1,2,3/_nfr_core.R          shared barrier ledger + pooled-ratio MDR estimator
code-semi/nofreeride-setting1,2,3/nfr_setting{1,2,3}_fig.R
code-semi/nofreeride-setting1,2,3/combine_nfr_settings_3x2.R
code-semi/nofreeride-setting1,2,3/nfr_compare_gam_dim.R
code-semi/nofreeride-setting1,2,3/nfr_samplesize_kde_ulsif_gam.R
code-semi/Figure2/Figure2_nfr.R
code-semi/sensity/nfr_D_n_vary.R
code-semi/setting1,2,3-appendix/nfr_setting*_varied_*.R
```

`code-semi/nofreeride` mirrors `nofreeride-setting1,2,3/`; several scripts source
`_nfr_core.R` through the shorter path.

Barrier scripts honour `NFR_REUSE_CACHE=1`, which reloads the saved `.rds` and
only redraws the figure:

```bash
NFR_REUSE_CACHE=1 Rscript "code-semi/nofreeride-setting1,2,3/nfr_setting1_fig.R"
```

The density-ratio ablation reads `NFR_REPS` and `NFR_NCORES` for optional preview
runs. Manuscript scripts otherwise keep their reported replication counts in the
source so a default run reproduces the stated experiment.

### Additional diagnostics (not figures in the manuscript)

```
code-semi/sensity/D_n_vary_setting2.R   the (n, beta, K0) sweep repeated under Setting (b)
code-semi/sensity/beta_error_diag.R     pi_t bias/MAE and Lmdr MAE as a function of beta
application/CCFD/CCFD_section5_revision.R   alpha-sweep summary table for the application
```

---

## Data

The application uses the **Credit Card Fraud Detection** dataset (284,807
transactions, 492 frauds; 28 PCA components plus `Time` and `Amount`).

1. Download `creditcard.csv` from
   <https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud>
2. Place it at `application/CCFD/creditcard.csv`
3. Verify: `md5sum` should be `e90efcb83d69faf99fcab8b0255024de`

`CCFD.R` falls back to `application/ccfd_data/creditcard.csv` if the file is not
found next to the script.

---

## Library reference — `OMDRC.R`

**Online decision rules**

| Function | Role |
|---|---|
| `OMDRC_OR(x.Lmdr, alpha, w, x.Lmdr_ini)` | Oracle Algorithm 1: consumes the true `Lmdr` and applies the local barrier. |
| `OMDRC_DD(z, z_ini, z0, z1, alpha, ..., w)` | Data-driven Algorithm 2: estimates the density ratio, tracks $\pi_t$ online, and applies the local barrier. |
| `OMDRC_OFF(x.Lmdr, alpha)` | Offline oracle fixed threshold (`FT` in the figures). |
| `OMDRC_FROM_RATIO(LR, LR_ini, alpha, ...)` | Shared engine, driven by precomputed density-ratio scores. |

All four return a list whose `de` element is the 0/1 decision vector.

**Density-ratio estimators** (all return an object consumed by `predict_ratio`)

| Function | Method |
|---|---|
| `fit_ratio_gam` | Probabilistic classifier via `mgcv` additive GAM — the estimator used in the paper (PC-DRE). |
| `fit_ratio_kde` | Kernel density ratio. |
| `fit_ratio_ulsif` | uLSIF direct density-ratio estimation. |
| `fit_ratio_classifier` | Penalised logistic classifier via `glmnet`. |

**Prevalence ($\pi_t$) estimators**

| Function | Role |
|---|---|
| `.estimate_local_pi_ratio` | Internal local mixture-MLE/EM update used by `OMDRC_DD`; its window includes the current observation. |
| `estimate_pi_online_ll` | Optional past-only local-linear diagnostic; not used by `OMDRC_DD`. |
| `estimate_pi_gam_batch` | Batch smoother. Offline benchmark only — not causal. |

**Baselines used in the comparisons**

| Function | Role |
|---|---|
| `STATIC_LMDR_DD` | Fixed threshold calibrated once on the warm-up batch. Loses MDR control when $\pi_t$ drifts. |
| `ROLLING_TOPK_DD` | Causal rolling top-$k$ policy. Its fixed $k$ can lose MDR control when prevalence rises. |

Adj-SAFFRON (from the `onlineFDR` package) is the online FDR comparator; `FT` is
the offline fixed-threshold oracle.

---

## Compute notes

* **Parallelism.** The simulation scripts use `foreach` + `doParallel` and claim
  all but two logical cores. Edit the `makeCluster(...)` call, or export
  `NFR_NCORES` for the barrier scripts, to change this.
* **Monte Carlo replications.** Main, sensitivity, and robustness experiments use
  the counts stated in their manuscript captions (normally 200; 150 for the
  density-ratio ablation and 100 for the high-dimensional DRE comparison).
* **Caching.** Each script saves its Monte Carlo output to a `*_results.rds`
  beside itself and reloads it when the corresponding reuse flag is on
  (`reuse_cache`/`REUSE_CACHE` in the script, or `NFR_REUSE_CACHE=1` in the
  environment). This makes figure/label tweaks cheap. The caches are not
  committed, so the first run of any script is a full run.
* **The application is the expensive one.** `CCFD.R` fits a 28-dimensional
  additive GAM density ratio and then sweeps $\alpha$; because result caches are
  not distributed, the first run recomputes and creates `ccfd_results_cache.rds`.

---

## Citation

```bibtex
@article{tian2026omdrc,
  title   = {Sequential Anomaly Detection with Online Missed Discovery Rate Control},
  author  = {Tian, Yang and Sun, Wenguang and Gang, Bowen},
  journal = {Statistics and Computing},
  year    = {2026}
}
```

## License

Released under the MIT License — see [LICENSE](LICENSE). The credit-card dataset
is **not** covered by this licence; it remains subject to its own terms on
Kaggle.
