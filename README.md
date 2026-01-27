# Sequential Anomaly Detection with Online Missed Discovery Rate Control

This repository contains the R implementation for the paper 'Sequential Anomaly Detection with Online Missed Discovery Rate Control' submitted to ICML 2026.

## 📖 Overview
OMDRC is a novel framework designed for high-stakes online monitoring where failing to detect a signal (Missed Discovery) is more costly than a false alarm. Unlike traditional Online FDR methods, our framework uses a dynamic "alpha-wealth" (earn-and-spend) mechanism to guarantee online MDR control, even in semi-supervised settings where signal patterns are learned on the fly.

---

## 🛠️ Prerequisites
The implementation is in **R**. Install the required packages using:

```R
install.packages(c("kedd", "onlineFDR", "ggplot2", "dplyr", "patchwork", "mvtnorm", "isoforest", "REBayes", "foreach", "doParallel"))
```

---

## 📊 Data Preparation
The real-world application in this paper uses the **ULB Credit Card Fraud dataset**.
1. Download `creditcard.csv` from [Kaggle Credit Card Fraud Detection](https://www.kaggle.com/datasets/mlg-ulb/creditcardfraud).
2. Place the file in a `data/` folder in the root directory (i.e., `./data/creditcard.csv`) before running the application script.

---

## 📂 Repository Structure

- `OMDRC.R`: **Core Algorithm Library**. Contains implementations for `OMDRC_OR` (Oracle), `OMDRC_DD` (Data-Driven), and other baselines.
- `setting1.R` / `setting2.R`: Main simulation scripts for Multi-modal and Skewed signal discovery.
- `capacity.R`: Analysis of the dynamic wealth allocation process under stationary vs. non-stationary (burst) regimes.
- `CCFD.R`: Script for the real-world credit card fraud detection application.
- `*_varied_*.R`: Scripts for sensitivity analyses of signal proportion ($\pi$), signal strength ($\mu, k$), and hyperparameters ($n, K_0, D, M$).

---

## 🚀 Reproducing Results

The following tables map the provided scripts to the figures in the paper.

### 1. Main Performance (MDR/FDR Control)
| Figure in Paper | Description | Script to Run |
| :--- | :--- | :--- |
| **Figure 1** | MDR Control vs. FDR Cost (Settings a & b) | `setting1_varied_pi.R` (with `pi_ = 0.08`) |
| **Figure 2** | Capacity Process & Signal Burst | `capacity.R` |
| **Figure 3** | Credit Card Fraud Detection | `CCFD.R` |

### 2. Sensitivity & Robustness (Appendix)
| Figure in Paper | Parameter Analyzed | Source Script |
| :--- | :--- | :--- |
| **Figure A.1 - A.4** | Signal Proportion ($\pi$) & Strength ($\mu, k$) | `setting1_varied_pi.R`, `setting2_varied_pi.R`, `setting1_varied_mu.R`, `setting2_varied_k.R` |
| **Figure A.5** | Hyperparameters ($n, K_0, D, M$) | `sample_sizes.R` |

---

## 💻 Usage Example

You can apply the Data-Driven OMDRC to your own data stream using the following logic:

```R
source("OMDRC.R")

# Parameters
alpha <- 0.1    # Target MDR level
D <- 1000       # Sliding window size for mixture density estimation
M <- 10         # Truncation threshold for the density ratio

# Run Data-Driven OMDRC
# z: online stream, z_ini: initial batch, z1: labeled alternative samples
results <- OMDRC_DD(z, z_ini, z1, alpha, D, M)

# Extract Decisions (1 = Reject/Signal, 0 = Not Reject/Null)
decisions <- results$de 
# Extract Density Ratio estimates
DR_estimates <- results$DR 
```

---

## 🔬 Core Mechanism: The $\alpha$-Wealth Process

The capacity (or "$\alpha$-wealth") $C_t$ evolves according to the following "earn-and-spend" logic using the estimated Density Ratio ($\widehat{\text{DR}}_t$):
- **Discovery (Reject):** $C_{t+1} = C_t + \alpha \cdot \widehat{\text{DR}}_t$
- **Non-Discovery:** $C_{t+1} = C_t - (1 - \alpha) \cdot \widehat{\text{DR}}_t$

This self-correcting mechanism ensures that the MDR is controlled at level $\alpha$ throughout the data stream, adapting to the evidence as it arrives.

---
