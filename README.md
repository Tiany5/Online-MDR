# Online Out-of-Distribution Testing with Missed Discovery Rate Control

This repository contains the official R implementation for the paper: **"Online Out-of-Distribution Testing with Missed Discovery Rate Control"**.

## 📖 Overview
OMDRC is a novel framework designed for high-stakes online monitoring where failing to detect a signal (Missed Discovery) is more costly than a false alarm. Unlike traditional Online FDR methods, our framework uses a dynamic "alpha-wealth" (earn-and-spend) mechanism to guarantees online MDR control.

---

## 🛠️ Prerequisites
The implementation is in **R**. Install the required packages using:

```R
install.packages(c("kedd", "onlineFDR", "ggplot2", "dplyr", "patchwork", "mvtnorm", "isoforest", "REBayes", "foreach", "doParallel"))
```

---

## 📂 Repository Structure

- `OMDRC.R`: **Core Algorithm Library**. Contains `OMDRC_OR` (Oracle), `OMDRC_DD` (Data-Driven), `OMDRC_OFF` (Offline).
- `setting1.R` / `setting2.R`: Main simulation scripts for Multi-modal and Skewed signal discovery with proportion $\pi=0.1$.
- `capacity.R`: Analysis of the dynamic wealth allocation process under stationary vs. non-stationary (burst) regimes.
- `*_vary.R`: Sensitivity analysis for $n$ (labeled size), $K_0$ (initial batch), $D$ (window size), $\pi$ (proportion), and $\mu$ (strength).

---

## 🚀 Reproducing Results

The following tables map the provided scripts to the figures in the paper.

### 1. Main Performance (MDR/FDR Control)
| Figure in Paper | Description | Script to Run |
| :--- | :--- | :--- |
| **Figure 1 (Left)** | MDR Control (Setting a & b) | `setting1.R`, `setting2.R` |
| **Figure 1 (Right)** | FDR Cost (Setting a & b) | `setting1.R`, `setting2.R` |
| **Figure 2** | Capacity Process & Signal Burst | `capacity.R` |
| **Figure 3** | Credit Card Fraud Detection | `CCFD.R` |

### 2. Sensitivity & Robustness (Appendix)
| Parameter Analyzed | Research Focus | Source Script |
| :--- | :--- | :--- |
| **Signal Proportion ($\pi$)** | Robustness to signal prevalence | `setting1_varied_pi.R`, `setting2_varied_pi.R` |
| **Signal Strength ($\mu, k$)** | Impact of signal-to-noise ratio | `setting1_varied_mu.R`, `setting2_varied_k.R` |
| **Hyperparameters ($n, K_0, D$)** | Influence of initialization and memory | `D,K0,n_vary.R` |

---

## 💻 Usage Example

You can apply the Data-Driven OMDRC to your own data stream using the following logic:

```R
source("OMDRC.R")

# Parameters
alpha <- 0.1    # Target MDR level
D <- 1000       # Sliding window size for null density estimation

# Run Data-Driven OMDRC
# z: online stream, z_ini: initial batch, z1: labeled alternative samples
results <- OMDRC_DD(z, z_ini, z1, alpha, D)

# Extract Decisions (1 = Reject/Signal, 0 = Not Reject/Null)
decisions <- results$de 
# Extract Density Ratios (DR estimates)
DR_estimates <- results$DR 
```

---

## 🔬 Core Mechanism: The $\alpha$-Wealth Process

The capacity $C_t$ evolves according to the following "earn-and-spend" logic:
- **Discovery (Reject):** $C_{t+1} = C_t + \alpha \cdot \text{Lmdr}_t$
- **Non-Discovery:** $C_{t+1} = C_t - (1 - \alpha) \cdot \text{Lmdr}_t$

This ensures that the empirical MDR is controlled at level $\alpha$ even when signal densities are estimated online.

---
