# ===================================================================
# 0. 基础设置
# ===================================================================
library(Matrix)
library(REBayes)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer) 
library(ggpubr)

setwd('/Users/tiany/Desktop/zju/online MDR/code-semi-github')
source('OMDRC.R')

# ===================================================================
# 1. 核心模拟函数 (修改为返回分子和分母矩阵)
# ===================================================================
run_simulation_ratio <- function(m, ini, n, pi, alpha, reps, D, methods_to_run = c("OR", "DD")) {
  
  # 初始化存储分子(Numerator)和分母(Denominator)的矩阵
  # 用于 MDR = E[Num]/E[Den]
  num.or <- matrix(0, reps, length(m)); den.or <- matrix(0, reps, length(m))
  num.dd <- matrix(0, reps, length(m)); den.dd <- matrix(0, reps, length(m))
  # FDR 仍按 E[V/R] 计算
  fdr.or <- matrix(0, reps, length(m)); fdr.dd <- matrix(0, reps, length(m))
  
  mu1 <- 2; mu2 <- -3; sd1 <- 0.7
  
  cl <- makeCluster(min(parallel::detectCores() - 1, 10))
  registerDoParallel(cl)
  clusterExport(cl, c("m","ini","n","pi", "reps", "alpha","D", "mu1","mu2","sd1", "methods_to_run"), envir = environment())
  
  result <- foreach(r = 1:reps, .packages = c("kedd","REBayes")) %dopar% {
    source('OMDRC.R')
    set.seed(r)
    
    N <- max(m) + ini
    theta <- rbinom(N, size = 1, prob = pi)
    z0 <- rnorm(N, 0, 1)
    mu_alt <- ifelse(rbinom(N, 1, 0.5) == 1, mu1, mu2)
    z_alt <- rnorm(N, mu_alt, sd1)
    z <- ifelse(theta == 0, z0, z_alt)
    
    f0 <- dnorm(z, 0, 1)
    f1 <- 0.5 * dnorm(z, mu1, sd1) + 0.5 * dnorm(z, mu2, sd1)
    z.lmdr <- (pi * f1) / ((1 - pi) * f0 + pi * f1)
    
    z_stream <- z[(ini + 1):N]
    z.lmdr_stream <- z.lmdr[(ini + 1):N]
    theta_stream <- theta[(ini + 1):N]
    
    # 生成标记样本
    z1 <- rnorm(n, ifelse(rbinom(n, 1, 0.5) == 1, mu1, mu2), sd1)
    
    res <- list()
    if ("OR" %in% methods_to_run) {
      de <- OMDRC_OR(z.lmdr_stream, alpha)$de
      res$num.or <- sapply(m, function(k) sum(theta_stream[1:k] * (1 - de[1:k])))
      res$den.or <- sapply(m, function(k) sum(theta_stream[1:k]))
      res$fdr.or <- sapply(m, function(k) sum((1 - theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
    }
    if ("DD" %in% methods_to_run) {
      de <- OMDRC_DD(z_stream, z[1:ini], z1, alpha, D)$de
      res$num.dd <- sapply(m, function(k) sum(theta_stream[1:k] * (1 - de[1:k])))
      res$den.dd <- sapply(m, function(k) sum(theta_stream[1:k]))
      res$fdr.dd <- sapply(m, function(k) sum((1 - theta_stream[1:k]) * de[1:k]) / max(sum(de[1:k]), 1))
    }
    res
  }
  stopCluster(cl)
  
  for (r in 1:reps) {
    if ("OR" %in% methods_to_run) {
      num.or[r,] <- result[[r]]$num.or; den.or[r,] <- result[[r]]$den.or; fdr.or[r,] <- result[[r]]$fdr.or
    }
    if ("DD" %in% methods_to_run) {
      num.dd[r,] <- result[[r]]$num.dd; den.dd[r,] <- result[[r]]$den.dd; fdr.dd[r,] <- result[[r]]$fdr.dd
    }
  }
  
  list(num.or=num.or, den.or=den.or, fdr.or=fdr.or, num.dd=num.dd, den.dd=den.dd, fdr.dd=fdr.dd)
}

# ===================================================================
# 2. 实验参数与执行
# ===================================================================
m_seq = seq(100, 1000, 50); pi_v = 0.1; alpha_v = 0.1; reps_v = 500; D_v = 1000
n_list <- c(5, 50, 100, 200, 500, 1000, 2000)
ini_list <- c(5, 50, 100, 200, 500, 1000, 2000)
D_list <- c(10, 50, 100, 200, 500, 1000, 2000)

# --- 基准: OMDRC.OR (计算一次) ---
base_res <- run_simulation_ratio(m_seq, 1000, 500, pi_v, alpha_v, reps_v, 1000, "OR")
mdr_or_vec <- colMeans(base_res$num.or) / pmax(colMeans(base_res$den.or), 1)
fdr_or_vec <- colMeans(base_res$fdr.or)

# --- 实验 A: 变化 n ---
res_n <- lapply(n_list, function(nv) run_simulation_ratio(m_seq, 1000, nv, pi_v, alpha_v, reps_v, 1000, "DD"))
# --- 实验 B: 变化 ini ---
res_ini <- lapply(ini_list, function(iv) run_simulation_ratio(m_seq, iv, 500, pi_v, alpha_v, reps_v, 1000, "DD"))
# --- 实验 C: 变化 D ---
res_D <- lapply(D_list, function(dv) run_simulation_ratio(m_seq, 1000, 500, pi_v, alpha_v, reps_v, dv, "DD"))

# ===================================================================
# 3. 数据处理与绘图函数
# ===================================================================
get_plot_df <- function(res_list, param_name, param_values, or_mdr, or_fdr) {
  df_or <- data.frame(t=m_seq, MDR=or_mdr, FDR=or_fdr, Method="OMDRC.OR")
  df_dd <- do.call(rbind, lapply(1:length(param_values), function(i) {
    res <- res_list[[i]]
    data.frame(t=m_seq, 
               MDR = colMeans(res$num.dd) / pmax(colMeans(res$den.dd), 1),
               FDR = colMeans(res$fdr.dd),
               Method = paste0("OMDRC.DD (", param_name, "=", param_values[i], ")"))
  }))
  df <- rbind(df_or, df_dd)
  df$Method <- factor(df$Method, levels = c("OMDRC.OR", unique(df_dd$Method)))
  return(df)
}

# --- 绘图配置 ---
custom_theme <- theme_bw() + theme(legend.position='bottom', legend.title=element_blank(), panel.border=element_rect(linewidth=0.8))
my_shapes <- c(16, 17, 15, 3, 4, 8, 2, 7)

# --- 绘图 A (n - Greens) ---
df_n <- get_plot_df(res_n, "n", n_list, mdr_or_vec, fdr_or_vec)
colors_n <- c("#F8766D", brewer.pal(9, "Greens")[3:9])
p_mdr_n <- ggplot(df_n, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_n) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(a.1)", y="MDR") + custom_theme
p_fdr_n <- ggplot(df_n, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_n) + scale_shape_manual(values=my_shapes) + labs(subtitle="(a.2)", y="FDR") + custom_theme

# --- 绘图 B (ini - Blues) ---
df_ini <- get_plot_df(res_ini, "K0", ini_list, mdr_or_vec, fdr_or_vec)
colors_ini <- c("#F8766D", brewer.pal(9, "Blues")[3:9])
p_mdr_ini <- ggplot(df_ini, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_ini) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(b.1)", y="MDR") + custom_theme
p_fdr_ini <- ggplot(df_ini, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_ini) + scale_shape_manual(values=my_shapes) + labs(subtitle="(b.2)", y="FDR") + custom_theme

# --- 绘图 C (D - Purples) ---
df_D <- get_plot_df(res_D, "D", D_list, mdr_or_vec, fdr_or_vec)
colors_D <- c("#F8766D", brewer.pal(9, "Purples")[3:9])
p_mdr_D <- ggplot(df_D, aes(t, MDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  geom_hline(yintercept=alpha_v, linetype="dashed") + scale_color_manual(values=colors_D) + 
  scale_shape_manual(values=my_shapes) + labs(subtitle="(c.1)", y="MDR") + custom_theme
p_fdr_D <- ggplot(df_D, aes(t, FDR, color=Method, shape=Method)) + geom_line() + geom_point() + 
  scale_color_manual(values=colors_D) + scale_shape_manual(values=my_shapes) + labs(subtitle="(c.2)", y="FDR") + custom_theme

# ===================================================================
# 4. 最终整合
# ===================================================================
fig_n <- (p_mdr_n | p_fdr_n) + plot_layout(guides="collect") & theme(legend.position="bottom")
fig_ini <- (p_mdr_ini | p_fdr_ini) + plot_layout(guides="collect") & theme(legend.position="bottom")
fig_D <- (p_mdr_D | p_fdr_D) + plot_layout(guides="collect") & theme(legend.position="bottom")

final_plot <- ggarrange(fig_n, fig_ini, fig_D, ncol=1, nrow=3)

print(final_plot)
# ggsave("MDR_Ratio_Expectations_Combined.pdf", final_plot, width=12, height=15)
