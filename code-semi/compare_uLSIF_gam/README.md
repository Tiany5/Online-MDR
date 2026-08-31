# compare_uLSIF — KDE vs uLSIF 密度比估计的样本量消融

本文件夹用于回应审稿意见：在**时变信号比例的在线 MDR 控制（OMDRC）**框架下，比较两种密度比 $r(x)=f_1(x)/f_0(x)$ 的估计方法——**KDE（核密度插件）** 与 **uLSIF（无约束最小二乘重要性拟合，直接密度比估计）**——并通过**参考样本量实验**证明二者的关系。

## 核心结论

> 本实验比较 KDE、uLSIF 与 GAM-based PC-DRE。PC-DRE 是正文 Algorithm 2
> 的默认估计器；KDE 在这些一维光滑设置中也能维持 MDR，而 uLSIF 在小、
> 中等参考样本量下偏差更明显。

机制：MDR 由容量规则对**逐点**局部 mdr 分数 $\hat q=\text{Lmdr}$ 的准确性决定；uLSIF 的目标函数在 $f_0$-加权 $L^2$ 范数下相合，对信号支撑集（$f_0$ 近零区域）的逐点精度收敛慢，故 MAE($\hat q$) 始终约为 KDE 的 3–4 倍。

## 文件说明

| 文件 | 说明 |
|------|------|
| `../OMDRC.R` | 共享的正文方法库，包含 PC-DRE、KDE、uLSIF、局部 $\pi_t$ 估计以及带 barrier 的 OMDRC。 |
| `sample_size_experiment.R` | 兼容旧文件名的入口；同样通过公共 `OMDRC_DD()` 使用 barrier。 |
| `../nofreeride-setting1,2,3/nfr_samplesize_kde_ulsif_gam.R` | 正文图的规范复现脚本。 |
| `../nofreeride-setting1,2,3/NFR_samplesize_kde_ulsif_gam_2x2_timeline.png/.pdf` | 正文结果图（2×2）。行 = setting，列 = MDR/FDR。 |

## 实验设置

两种代表性数据几何，各配一条 $\pi_t$ 路径：

| 图行 | Setting | $F_0$ / $F_1$ | $\pi_t$ 路径 |
|------|---------|---------------|-------------|
| (a) | Multimodal | $N(0,1)$ / $0.5N(-3,0.7^2)+0.5N(2,0.7^2)$ | 单调下降（prior drop）0.30→0.08 |
| (b) | Gaussian shift | $N(0,1)$ / $N(2.5,1)$ | 非单调突起（surge）0.05→峰0.30→0.05 |

图列：`.1` = MDR，`.2` = FDR，即 (a.1)(a.2)(b.1)(b.2)。

## 关键参数

- 参考样本量：$n_0=n_1\in\{500,1000,2000,4000\}$
- 在线流 $T=1000$；`K0=500`，`D_min=10`，`beta=0.6`，`w=100`，`alpha=0.1`，`reps=150`（与图注一致）
- **uLSIF 容量随样本量增长**（满足非参一致性）：中心数 $\propto\sqrt{n}$（`n_centers_schedule`），带宽网格向更细延伸，岭参数向更弱延伸
- `gs_mean=2.5`：Gaussian shift 的信号强度（$F_1=N(2.5,1)$），使 uLSIF 恰好在 n≈2000 处跨过 $\alpha$

规范脚本会同时输出逐时 MDR/FDR、95% Monte Carlo 置信带，以及终端
MDR、MAE($\widehat q$) 和 RMSE($\widehat\pi$) 汇总，避免 README 中的静态
数值与重新运行的结果缓存脱节。

## 可视化规范

- KDE = 橘色系，uLSIF = 蓝色系；每族内按 $n$ 由浅到深；点形状随 $n$（●▲■◆）；color 与 shape 映射同一分组以合并为单一图例。
- $\alpha=0.1$ 以黑色虚线标注于 MDR 面板。

## 运行方式

```bash
Rscript "code-semi/nofreeride-setting1,2,3/nfr_samplesize_kde_ulsif_gam.R"
```

依赖：`foreach`、`doParallel`、`ggplot2`、`dplyr`、`tidyr`、`patchwork`。

### 缓存复用（按需重算）

默认会重算并写入带 `NFR_` 前缀的缓存。仅重绘时可设置
`NFR_REUSE_CACHE=1`；预览运行可通过 `NFR_REPS` 和 `NFR_NCORES` 降低重复数
和并行进程数，预览输出会使用不同文件名，不会覆盖正文结果。

## 备注 / 边界

uLSIF 等直接密度比估计的真正优势在**高维**场景（此时 KDE 受维度灾难影响）。本文结论限定于**低维光滑**的流式检验设置，在此范围内 KDE 更高效。
