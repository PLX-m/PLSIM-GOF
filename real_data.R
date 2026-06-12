# Spline single-index
library(splines2)
library(MASS)
library(openxlsx)
library(caret)
library(glmnet)
library(randomForest)
library(reticulate)

setwd("code/")
source("SIM_functions_final.R")

# 指定 Python 环境
use_python("D:/Anaconda_envs/envs/myenv/python.exe", required=TRUE)
# 加载 Python 文件
source_python("code/python_code.py")



PLSIM_test_split_ratio_cross2 <- function(
    Y, X, Z,
    split_ratio = c(4, 4),
    K_si = 4,
    C_lambda = 0.5,
    mu_true = NULL,
    beta_true = NULL,
    oracle_si = NULL,
    seed = 2026
){
  
  set.seed(seed)
  
  Z <- as.matrix(Z)
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  
  n <- nrow(Z)
  p <- ncol(X)
  q <- ncol(Z)
  
  eps_num <- 1e-12
  
  stopifnot(length(Y) == nrow(X), nrow(X) == nrow(Z))
  stopifnot(length(split_ratio) == 2)
  stopifnot(all(split_ratio > 0))
  
  np <- reticulate::import("numpy")
  
  ############################################################
  ## 0. 先把全部数据等分成两份
  ############################################################
  
  id_all <- sample(seq_len(n))
  
  n_1 <- floor(n / 2)
  
  id_part1 <- id_all[seq_len(n_1)]
  id_part2 <- id_all[(n_1 + 1):n]
  
  stopifnot(length(id_part1) > 20, length(id_part2) > 20)
  
  
  ############################################################
  ## 1. 单个 half 内部完整走一遍流程：
  ##
  ##   1) 按 split_ratio 切成 block1 和 id_rf
  ##   2) block1 再切成 id_ce 和 id_beta
  ##   3) id_ce 估计条件期望
  ##   4) id_beta 估计 beta
  ##   5) id_rf 内部 K 折 OOF 估计单指标 residual
  ##   6) 用 id_rf 的 OOF residual 训练 RF
  ############################################################
  
  run_one_half <- function(id_half, seed) {
    
    id_half <- sample(id_half)  # 随机打乱这一半数据里的样本编号顺序
    
    n_half <- length(id_half)  # 221
    
    ratio_use <- split_ratio / sum(split_ratio)  # 把传入的比例转换成总和为 1 的比例。
    
    # 第一大块：ce + beta
    # 第二大块：rf
    n_block1 <- floor(n_half * ratio_use[1])  # 110
    n_rf     <- n_half - n_block1  # 111
    
    # 第一大块再对半分
    n_ce   <- floor(n_block1 / 2)
    n_beta <- n_block1 - n_ce
    
    stopifnot(n_ce > 5, n_beta > 5, n_rf > 10)
    
    group_id <- rep(
      c("ce", "beta", "rf"),
      times = c(n_ce, n_beta, n_rf)
    )
    
    id_split <- split(id_half, group_id)
    
    id_ce   <- id_split$ce
    id_beta <- id_split$beta
    id_rf   <- id_split$rf
    
    
    ############################################################
    ## 1.1 用 id_ce 估计条件期望；id_beta 估计 beta
    ############################################################
    
    if (!is.null(beta_true)) {
      
      beta_use <- as.numeric(beta_true)
      
    } else {
      
      if (!is.null(mu_true)) {
        
        E_YZ_beta <- as.numeric(mu_true$mu_Y)
        E_XZ_beta <- mu_true$mu_X
        
        Y_tilde_beta <- Y[id_beta] - E_YZ_beta[id_beta]
        
        X_tilde_beta <- sweep(
          X[id_beta, , drop = FALSE], 1, E_XZ_beta[id_beta], "-")
        
      } else {
        
        fit_ce <- do.call(
          fit_conditional_expectation_models,
          list(
            Z_train = np$array(Z[id_ce, , drop = FALSE]),
            Y_train = np$array(Y[id_ce]),
            X_train = np$array(X[id_ce, , drop = FALSE]),
            verbose = FALSE,
            seed = as.integer(seed + 100)
          )
        )
        
        res_beta <- do.call(
          predict_conditional_expectations,
          list(
            fit_obj = fit_ce,
            Z_new = np$array(Z[id_beta, , drop = FALSE])
          )
        )
        
        E_YZ_beta <- as.numeric(res_beta[[1]])
        E_XZ_beta <- as.matrix(res_beta[[2]])
        
        Y_tilde_beta <- Y[id_beta] - E_YZ_beta
        X_tilde_beta <- X[id_beta, , drop = FALSE] - E_XZ_beta
      }
      
      
      n_beta_obs <- length(Y_tilde_beta)
      d_beta     <- ncol(X_tilde_beta)
      
      if (d_beta > n_beta_obs) {
        
        # 高维情形：d > n，用 Lasso
        nfolds_beta <- min(5, n_beta_obs)
        
        fit_cv <- glmnet::cv.glmnet(
          x = X_tilde_beta,
          y = Y_tilde_beta,
          alpha = 1,
          intercept = FALSE,
          standardize = TRUE,
          nfolds = nfolds_beta
        )
        
        beta_use <- as.numeric(coef(fit_cv, s = "lambda.min")[-1])
        
      } else {
        
        # 低维情形：d <= n，用 OLS
        fit_ols <- lm.fit(
          x = X_tilde_beta,
          y = Y_tilde_beta
        )
        
        beta_use <- as.numeric(fit_ols$coefficients)
        beta_use[is.na(beta_use)] <- 0
      }
      
    }
    
    
    ############################################################
    ## 1.2 在 id_rf 内部做 OOF 单指标 residual
    ############################################################
    
    n_rf_actual <- length(id_rf)
    K_si_use <- min(K_si, n_rf_actual)
    
    fold_si <- sample(rep(seq_len(K_si_use), length.out = n_rf_actual))
    
    r_rf_oof <- rep(NA_real_, n_rf_actual)
    
    for (k in seq_len(K_si_use)) {
      
      id_si_test_local  <- which(fold_si == k)
      id_si_train_local <- which(fold_si != k)
      
      id_si_test_global  <- id_rf[id_si_test_local]
      id_si_train_global <- id_rf[id_si_train_local]
      
      if (!is.null(oracle_si)) {
        
        Rhat_si_test <- oracle_si$g_fun(
          as.vector(
            Z[id_si_test_global, , drop = FALSE] %*% oracle_si$theta
          )
        )
        
      } else {
        
        R_si_train <- Y[id_si_train_global] -
          as.vector(
            X[id_si_train_global, , drop = FALSE] %*% beta_use
          )
        
        fit_si_k <- spline_single_index(
          Z_train = Z[id_si_train_global, , drop = FALSE],
          Y = R_si_train
        )
        
        Rhat_si_test <- predict_merged_si(
          theta_merged = fit_si_k$theta_hat_orig,
          g_coef_merged = fit_si_k$g_coef_hat,
          Z_new = Z[id_si_test_global, , drop = FALSE],
          Z_mu = fit_si_k$Z_mu,
          Z_sd = fit_si_k$Z_sd,
          N = fit_si_k$N,
          a = fit_si_k$a,
          d = ncol(Z)
        )
      }
      
      r_rf_oof[id_si_test_local] <-
        Y[id_si_test_global] -
        as.vector(
          X[id_si_test_global, , drop = FALSE] %*% beta_use) - Rhat_si_test
    }
    
    valid_rf <- is.finite(r_rf_oof)
    
    id_rf_valid <- id_rf[valid_rf]
    
    X_rf <- X[id_rf_valid, , drop = FALSE]
    Z_rf <- Z[id_rf_valid, , drop = FALSE]
    r_rf <- r_rf_oof[valid_rf]
    
    
    ############################################################
    ## 1.3 用这个 half 的 OOF residual 训练 RF
    ############################################################
    
    X_mu <- colMeans(X_rf)
    X_sd <- apply(X_rf, 2, sd)
    X_sd[X_sd == 0 | !is.finite(X_sd)] <- 1
    
    Z_mu <- colMeans(Z_rf)
    Z_sd <- apply(Z_rf, 2, sd)
    Z_sd[Z_sd == 0 | !is.finite(Z_sd)] <- 1
    
    X_rf_std <- sweep(sweep(X_rf, 2, X_mu, "-"), 2, X_sd, "/")
    Z_rf_std <- sweep(sweep(Z_rf, 2, Z_mu, "-"), 2, Z_sd, "/")
    
    rf_final <- randomForest::randomForest(
      x = as.data.frame(cbind(X_rf_std, Z_rf_std)),
      y = r_rf,
      ntree = 500,
      mtry = max(5, floor(sqrt(ncol(X_rf_std) + ncol(Z_rf_std)))),
      nodesize = 5
    )
    
    
    ############################################################
    ## 1.4 返回该 half 的完整结果
    ############################################################
    
    list(
      beta_use = beta_use,
      
      id_ce = id_ce,
      id_beta = id_beta,
      id_rf = id_rf,
      id_rf_valid = id_rf_valid,  # 训练RF的数据编号
      
      r_rf = r_rf,          # OOF 残差
      rf_final = rf_final,  # RF 模型
      
      X_mu = X_mu,          # 训练RF用到的X均值
      X_sd = X_sd,          # 训练RF用到的X标准差
      Z_mu = Z_mu,
      Z_sd = Z_sd,
      
      n_ce = length(id_ce),
      n_beta = length(id_beta),
      n_rf = length(id_rf),
      n_rf_valid = length(id_rf_valid),
      K_si = K_si_use
    )
  }
  
  
  ############################################################
  ## 2. 两份数据都完整走一遍
  ############################################################
  
  fit_part1 <- run_one_half(
    id_half = id_part1, seed = seed
  )
  
  fit_part2 <- run_one_half(
    id_half = id_part2, seed = seed
  )
  
  
  ############################################################
  ## 3. 用 one half 的 RF 去预测另一个 half 的 id_rf_valid，
  ##    然后和另一个 half 自己得到的 OOF residual 算 Zn 和 p 值
  ############################################################
  
  calc_Zn_between_halves <- function(train_fit, test_fit) {
    
    id_test <- test_fit$id_rf_valid
    
    X_test <- X[id_test, , drop = FALSE]
    Z_test <- Z[id_test, , drop = FALSE]
    
    r_test <- as.numeric(test_fit$r_rf)
    
    ############################################################
    ## 1. 用 train_fit 的标准化参数，把 test half 标准化
    ##    注意：这里必须用 train_fit 的 X_mu, X_sd, Z_mu, Z_sd，
    ##    因为 RF 是在 train_fit 的标准化尺度上训练的。
    ############################################################
    
    X_test_std <- sweep(
      sweep(X_test, 2, train_fit$X_mu, "-"), 2, train_fit$X_sd, "/")
    
    Z_test_std <- sweep(
      sweep(Z_test, 2, train_fit$Z_mu, "-"), 2, train_fit$Z_sd, "/")
    
    ############################################################
    ## 2. RF 预测另一个 half 的 residual predictability
    ############################################################
    
    f_pred <- as.numeric(
      predict(
        train_fit$rf_final,
        newdata = as.data.frame(cbind(X_test_std, Z_test_std))
      )
    )
    
    valid <- is.finite(r_test) & is.finite(f_pred) &
      apply(X_test_std, 1, function(v) all(is.finite(v)))
    
    r_use <- r_test[valid]
    f_use <- f_pred[valid]
    X_use <- X_test_std[valid, , drop = FALSE]
    
    n_use <- length(r_use)
    p_use <- ncol(X_use)
    
    if (n_use <= 5) {
      return(list(
        Zn = NA_real_,
        p_value = NA_real_,
        n_use = n_use,
        f_pred = f_pred,
        f_orth = rep(NA_real_, length(f_pred)),
        r_test = r_test,
        gamma_hat = rep(NA_real_, p_use),
        Orthogonal_beta = NA_integer_
      ))
    }
    
    
    ############################################################
    ## 3. 线性正交化：
    ##    用 Lasso 把 f_pred 中可由 X 线性解释的部分去掉
    ##
    ##    f_orth = f_pred - X gamma_hat
    ##
    ##    用 square-root-lasso 选 lambda：
    ##    lambda_theory = C_lambda * sqrt(log(p) / n)
    ##    再乘上每个 lambda 下 residual 的尺度 scale_vec
    ############################################################
    
    fit_gam <- glmnet::glmnet(
      x = X_use,
      y = f_use,
      alpha = 1,
      intercept = FALSE,
      standardize = FALSE,
      nlambda = 100
    )
    

    W_mat <- f_use - X_use %*% as.matrix(fit_gam$beta)
    
    scale_vec <- apply(W_mat, 2, function(w) {
      sqrt(sum(w^2)) / sqrt(length(w))
    })
    
    lambda_theory <- C_lambda * sqrt(log(p_use) / n_use)
    
    index <- which.min(abs(fit_gam$lambda - lambda_theory * scale_vec))
    
    gamma_hat <- as.numeric(fit_gam$beta[, index])
    
    Orthogonal_beta <- sum(abs(gamma_hat) > 1e-8)
    
    f_orth_use <- f_use - as.vector(X_use %*% gamma_hat)
    
    
    ############################################################
    ## 4. 用正交化后的 f_orth 和 residual r 计算统计量
    ############################################################
    
    r_c <- r_use - mean(r_use)
    f_c <- f_orth_use - mean(f_orth_use)
    
    rw <- r_c * f_c
    S <- sum(rw) / sqrt(n_use)
    V <- mean(rw^2)
    Zn <- S / sqrt(V + eps_num)
    p_value <- 1 - pnorm(Zn)

  
    f_orth <- rep(NA_real_, length(f_pred))
    f_orth[valid] <- f_orth_use
    
    list(
      Zn = Zn,
      p_value = p_value,
      n_use = n_use,
      
      f_pred = f_pred,
      f_orth = f_orth,
      r_test = r_test,
      
      gamma_hat = gamma_hat,
      Orthogonal_beta = Orthogonal_beta,
      lambda_theory = lambda_theory,
      lambda_selected = fit_gam$lambda[index],
      scale_selected = scale_vec[index]
    )
  }
  
  ############################################################
  ## 4. 两个方向：
  ##    part1 的 f 预测 part2 的 r
  ##    part2 的 f 预测 part1 的 r
  ############################################################
  
  stat_12 <- calc_Zn_between_halves(
    train_fit = fit_part1,
    test_fit = fit_part2
  )
  
  stat_21 <- calc_Zn_between_halves(
    train_fit = fit_part2,
    test_fit = fit_part1
  )
  
  
  ############################################################
  ## 5. CCT 组合两个 p 值
  ############################################################
  
  cct_pvalue <- function(pvals, weights = NULL) {
    
    pvals <- as.numeric(pvals)
    pvals <- pvals[is.finite(pvals)]
    
    if (length(pvals) == 0) {
      return(NA_real_)
    }
    
    if (is.null(weights)) {
      weights <- rep(1 / length(pvals), length(pvals))
    } else {
      weights <- weights / sum(weights)
    }
    
    pvals <- pmin(pmax(pvals, 1e-15), 1 - 1e-15)
    
    T_cct <- sum(weights * tan((0.5 - pvals) * pi))
    
    p_cct <- 0.5 - atan(T_cct) / pi
    
    pmin(pmax(p_cct, 0), 1)
  }
  
  p_CCT <- cct_pvalue(
    pvals = c(stat_12$p_value, stat_21$p_value)
  )
  
  
  ############################################################
  ## 6. 输出
  ############################################################
  
  list(
    Zn_12 = stat_12$Zn,
    p_12 = stat_12$p_value,
    n_12 = stat_12$n_use,
    Orthogonal_beta_12 = stat_12$Orthogonal_beta,
    lambda_selected_12 = stat_12$lambda_selected,
    
    Zn_21 = stat_21$Zn,
    p_21 = stat_21$p_value,
    n_21 = stat_21$n_use,
    Orthogonal_beta_21 = stat_21$Orthogonal_beta,
    lambda_selected_21 = stat_21$lambda_selected,
    
    p_CCT = p_CCT,
    
    id_part1 = id_part1,
    id_part2 = id_part2,
    
    fit_part1 = fit_part1,
    fit_part2 = fit_part2,
    
    stat_12 = stat_12,
    stat_21 = stat_21
  )
}



##################################################################################################
start_time <- Sys.time()


library(lars)
data(diabetes)
dim(diabetes)   # 442 3

X <- diabetes$x   # 442 x 10
Y <- diabetes$y   # length 442

dim(X)        # 442 10
length(Y)     # 442
colnames(X)   # 已经有了：age sex bmi map tc ldl hdl tch ltg glu

X0 <- diabetes$x
Y  <- diabetes$y

# 442   2
Z <- X0[, c("bmi", "age", "map"), drop = FALSE]  # bp 在这个对象里叫 map
# 442   8
X <- X0[, setdiff(colnames(X0), c("bmi", "age","map")), drop = FALSE]

# 变量标准化
X <- scale(X, center = TRUE, scale = TRUE)
Z <- scale(Z, center = TRUE, scale = TRUE)


result = PLSIM_test_split_ratio_cross2(
  Y=Y, X=X, Z=Z, 
  split_ratio = c(4, 4),
  K_si = 4,
  C_lambda = 0.5
)


result$Zn_12   # -0.6951736
result$p_12    # 0.7565268

result$Zn_21   # -1.296672
result$p_21    # 0.902628

result$p_CCT   # 0.85878


end_time = Sys.time()
end_time - start_time  # 27.17573 secs

# 与2013年以及2024年文章得出的结论相同，都认为数据服从PLSIM。
# 但是2013是在X和Y被变量U扭曲情况下得出的结论，我们直接用X和Z，没用U扭曲数据。
# 2014年像我们一样，不考虑变量U扭曲，他直接把bmi放到单指标变量里面。

