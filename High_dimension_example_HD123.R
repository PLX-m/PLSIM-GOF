# Spline single-index
library(splines2)
library(MASS)
library(openxlsx)
library(caret)
library(glmnet)
library(randomForest)
library(reticulate)

source("code/SIM_functions.R")
source("code/Orthogonalization.R")

# Specify the Python interpreter; change this path to your local Python environment.
use_python("D:/Anaconda_envs/envs/myenv/python.exe", required = TRUE)
# 加载 Python 文件
source_python("code/python_code.py")


PLSIM_test_split_ratio_cross2 <- function(
    Y, X, Z, 
    split_ratio = c(4, 4),
    K_si = 4,
    C_lambda_theta = 0.3,
    C_lambda_beta = 1,
    df_si_orth = 6,
    do_si_orth = TRUE,
    beta_true = NULL,
    oracle_si = NULL,
    seed = 1
) {
  
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
  
  if (K_si < 2L) stop("K_si must be at least 2.")
  
  np <- reticulate::import("numpy")
  
  
  ############################################################
  ## 0. 先把全部数据等分成两份，并各自划分
  ############################################################
  
  id_all <- sample(seq_len(n))
  n_1 <- floor(n / 2)
  
  id_part1 <- id_all[seq_len(n_1)]
  id_part2 <- id_all[(n_1 + 1):n]
  
  stopifnot(length(id_part1) > 20, length(id_part2) > 20)
  
  
  ############################################################
  ## 1. 先在两个 half 内部分别划分 ce / beta / rf
  ############################################################
  
  make_half_split <- function(id_half) {
    
    id_half <- sample(id_half)
    n_half <- length(id_half)
    ratio_use <- split_ratio / sum(split_ratio)
    
    n_block1 <- floor(n_half * ratio_use[1])
    n_rf <- n_half - n_block1
    
    n_ce <- floor(n_block1 / 2)
    n_beta <- n_block1 - n_ce
    
    stopifnot(n_ce > 5, n_beta > 5, n_rf > 10)
    
    group_id <- rep(c("ce", "beta", "rf"), times = c(n_ce, n_beta, n_rf))
    id_split <- split(id_half, group_id)
    
    list(
      id_ce = id_split$ce,
      id_beta = id_split$beta,
      id_rf = id_split$rf
    )
  }
  
  split1 <- make_half_split(id_part1)
  split2 <- make_half_split(id_part2)
  
  
  id_si_scale <- c(
    split1$id_ce, split1$id_beta,
    split2$id_ce, split2$id_beta
  )
  
  Z_scale_ref <- Z[id_si_scale, , drop = FALSE]
  Z_mu_si <- colMeans(Z_scale_ref)
  Z_sd_si <- apply(Z_scale_ref, 2, sd)
  
  if (any(!is.finite(Z_mu_si)) || any(!is.finite(Z_sd_si)) || any(Z_sd_si <= 0)) {
    stop("Invalid common Z scaler for single-index fitting.")
  }
  
  Z_scale_std <- sweep(sweep(Z_scale_ref, 2, Z_mu_si, "-"), 2, Z_sd_si, "/")
  a_si <- as.numeric(
    quantile(sqrt(rowSums(Z_scale_std^2)), probs = 0.95, names = FALSE)
  )
  
  
  ############################################################
  ## 2. 单个 half 内部完整走一遍流程
  ############################################################
  
  run_one_half <- function(split_obj, seed) {
    
    id_ce <- split_obj$id_ce
    id_beta <- split_obj$id_beta
    id_rf <- split_obj$id_rf
    
    ############################################################
    ## 1.1 用 id_ce 估计条件期望；id_beta 估计 beta
    ############################################################
    
    if (!is.null(beta_true)) {
      
      beta_use <- as.numeric(beta_true)
      
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
      
      
      n_beta_obs <- length(Y_tilde_beta)
      d_beta <- ncol(X_tilde_beta)
      
      if (d_beta > n_beta_obs) {
        
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
        
        fit_ols <- lm.fit(x = X_tilde_beta, y = Y_tilde_beta)
        
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
    
    N_si_fold <- rep(NA_integer_, K_si_use)   
    N_si_half <- NA_integer_                 
    si_fit_list <- vector("list", K_si_use)
    
    for (k in seq_len(K_si_use)) {
      
      id_si_test_local <- which(fold_si == k)
      id_si_train_local <- which(fold_si != k)
      
      id_si_test_global <- id_rf[id_si_test_local]
      id_si_train_global <- id_rf[id_si_train_local]
      
      if (!is.null(oracle_si)) {
        
        # 计算 oracle fitted value
        index_si_test <- as.vector(
          Z[id_si_test_global, , drop = FALSE] %*% oracle_si$theta
        )
        
        Rhat_si_test <- oracle_si$g_fun(index_si_test)
        
      } else {
        
        R_si_train <- Y[id_si_train_global] - 
          as.vector(X[id_si_train_global, , drop = FALSE] %*% beta_use)
        
        fit_si_k <- spline_single_index(
          Z_train = Z[id_si_train_global, , drop = FALSE],
          Y = R_si_train,
          Z_mu = Z_mu_si,
          Z_sd = Z_sd_si,
          a_si = a_si
        )
        
        # 后面在 opposite evaluation half 的 Z 上计算 u 和 S
        si_fit_list[[k]] <- fit_si_k
        
        N_si_fold[k] <- as.integer(fit_si_k$N)   # same
        if (is.na(N_si_half)) {
          N_si_half <- as.integer(fit_si_k$N)
        }
        
        Rhat_si_test <- predict_si(
          theta = fit_si_k$theta_hat,
          g_coef = fit_si_k$g_coef_hat,
          Z_new = Z[id_si_test_global, , drop = FALSE],
          Z_mu = fit_si_k$Z_mu,
          Z_sd = fit_si_k$Z_sd,
          N = fit_si_k$N,
          a = fit_si_k$a,
          d = ncol(Z)
        )
      }
      
      r_rf_oof[id_si_test_local] <- Y[id_si_test_global] - as.vector(
        X[id_si_test_global, , drop = FALSE] %*% beta_use) - Rhat_si_test
    }
    
    
    if (is.null(oracle_si)) {
      
      N_si_unique <- unique(N_si_fold[is.finite(N_si_fold)])
      
      if (length(N_si_unique) > 1L) {
        warning(
          "Different spline dimensions were selected across OOF folds. ",
          "The pooled SI orthogonalization uses N_si = ",
          N_si_half,
          "."
        )
      }
    }
    
    if (any(!is.finite(r_rf_oof))) {
      stop(sprintf(
        "Non-finite OOF residuals detected: %d.",
        sum(!is.finite(r_rf_oof))
      ))
    }
    
    id_rf_valid <- id_rf
    
    X_rf <- X[id_rf_valid, , drop = FALSE]
    Z_rf <- Z[id_rf_valid, , drop = FALSE]
    r_rf <- r_rf_oof
    fold_si_rf <- fold_si
    
    
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
      mtry = max(1, floor(sqrt(ncol(X_rf_std) + ncol(Z_rf_std)))),
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
      id_rf_valid = id_rf_valid,
      
      r_rf = r_rf,
      rf_final = rf_final,
      
      fold_si_rf = fold_si_rf,
      si_fit_list = si_fit_list,
      Z_mu_si = Z_mu_si,
      Z_sd_si = Z_sd_si,
      
      X_mu = X_mu,
      X_sd = X_sd,
      Z_mu = Z_mu,
      Z_sd = Z_sd,
      
      n_ce = length(id_ce),
      n_beta = length(id_beta),
      n_rf = length(id_rf),
      n_rf_valid = length(id_rf_valid),
      K_si = K_si_use,
      N_si = N_si_half,
      N_si_fold = N_si_fold
    )
  }
  
  
  ############################################################
  ## 2. 两份数据都完整走一遍
  ############################################################
  
  fit_part1 <- run_one_half(split_obj = split1, seed = seed)
  fit_part2 <- run_one_half(split_obj = split2, seed = seed)
  
  
  ############################################################
  ## 3. 一个 half 的 RF 预测另一个 half
  ############################################################
  
  calc_Zn_between_halves <- function(train_fit, test_fit) {
    
    id_test <- test_fit$id_rf_valid
    
    X_test <- X[id_test, , drop = FALSE]
    Z_test <- Z[id_test, , drop = FALSE]
    
    r_test <- as.numeric(test_fit$r_rf)
    fold_test <- as.integer(test_fit$fold_si_rf)
    
    # 正交化方向必须由 opposite half 的第 k 个 SI fit 在当前 test half 的第 k 折协变量上重新评价
    if (do_si_orth && is.null(oracle_si)) {
      
      if (train_fit$K_si != test_fit$K_si) {
        stop("The two half-samples must use the same number of SI folds.")
      }
      
      u_test <- rep(NA_real_, length(id_test))
      theta_score_test <- matrix(
        NA_real_,
        nrow = length(id_test),
        ncol = q
      )
      
      for (k in seq_len(test_fit$K_si)) {
        
        idx_k <- which(fold_test == k)
        if (length(idx_k) == 0L) next
        
        fit_si_ref <- train_fit$si_fit_list[[k]]
        if (is.null(fit_si_ref)) {
          stop(sprintf("Missing opposite-half SI fit for fold %d.", k))
        }
        
        score_k <- build_oof_theta_score(
          Z_new = Z_test[idx_k, , drop = FALSE],
          fit_si = fit_si_ref,
          eps = 1e-12
        )
        
        u_test[idx_k] <- score_k$u
        theta_score_test[idx_k, ] <- score_k$theta_score
      }
      
    } else {
      
      # oracle / 不做 SI 正交化时，仅用于保持后续长度一致
      u_test <- rep(0, length(id_test))
      theta_score_test <- matrix(0, nrow = length(id_test), ncol = q)
    }
    
    ############################################################
    ## 1. 用 train_fit 的参数标准化 test half
    ############################################################
    
    X_test_std <- sweep(
      sweep(X_test, 2, train_fit$X_mu, "-"), 2, train_fit$X_sd, "/"
    )
    
    Z_test_std <- sweep(
      sweep(Z_test, 2, train_fit$Z_mu, "-"), 2, train_fit$Z_sd, "/"
    )
    
    
    ############################################################
    ## 2. RF 预测 residual predictability
    ############################################################
    
    f_pred <- as.numeric(
      predict(
        train_fit$rf_final,
        newdata = as.data.frame(cbind(X_test_std, Z_test_std))
      )
    )
    
    if (any(!is.finite(r_test))) {
      stop("Non-finite test residual.")
    }
    
    valid <- is.finite(f_pred) & is.finite(u_test) &
      apply(X_test_std, 1, function(v) all(is.finite(v)))
    
    if (do_si_orth && is.null(oracle_si)) {
      valid <- valid & apply(theta_score_test, 1, function(v) all(is.finite(v)))
    }
    
    r_use <- r_test[valid]
    f_use <- f_pred[valid]
    u_use <- u_test[valid]
    
    X_use <- X_test_std[valid, , drop = FALSE]
    theta_score_use <- theta_score_test[valid, , drop = FALSE]
    fold_use <- fold_test[valid]
    
    n_use <- length(r_use)
    p_use <- ncol(X_use)
    
    if (n_use <= 5) {
      stop(sprintf("Effective sample size is too small: n_use = %d.", n_use))
    }
    
    
    ############################################################
    ## 3. 线性正交化辅助函数
    ############################################################
    
    linear_orthogonalize <- function(f_input) {
      
      fit_gam <- glmnet::glmnet(
        x = X_use,
        y = f_input,
        alpha = 1,
        intercept = TRUE,
        standardize = FALSE,
        nlambda = 100
      )
      
      W_mat <- as.matrix(X_use %*% fit_gam$beta)
      W_mat <- sweep(W_mat, 2, fit_gam$a0, "+")
      W_mat <- sweep(W_mat, 1, f_input, function(pred, obs) obs - pred)
      
      scale_vec <- apply(W_mat, 2, function(w) sqrt(mean(w^2)))
      scale_vec[!is.finite(scale_vec) | scale_vec < eps_num] <- eps_num
      
      lambda_theory <- C_lambda_beta * sqrt(log(max(p_use, 2L)) / n_use)
      
      index <- which.min(
        abs(fit_gam$lambda - lambda_theory * scale_vec)
      )
      
      gamma_hat <- as.numeric(fit_gam$beta[, index])
      alpha_hat <- as.numeric(fit_gam$a0[index])
      f_orth <- f_input - alpha_hat - as.vector(X_use %*% gamma_hat)
      
      linear_orth_before <- max(
        abs(as.vector(crossprod(X_use, f_input)) / n_use)
      )
      
      linear_orth_after <- max(
        abs(as.vector(crossprod(X_use, f_orth)) / n_use)
      )
      
      list(
        f_orth = f_orth,
        linear_orth_before = linear_orth_before,
        linear_orth_after = linear_orth_after
      )
    }
    
    ############################################################
    ## 4. 最终方法：spline + linear + theta-score 联合正交化
    ############################################################
    
    spline_orth_before <- NA_real_
    spline_orth_after <- NA_real_
    theta_orth_before <- NA_real_
    theta_orth_after <- NA_real_
    linear_orth_before_final <- NA_real_
    linear_orth_after_final <- NA_real_
    
    path_index <- NA_real_
    path_ratio <- NA_real_
    path_gap <- NA_real_
    linear_kkt_ratio <- NA_real_
    theta_kkt_ratio <- NA_real_
    theta_kkt_ratio_by_fold <- rep(NA_real_, test_fit$K_si)
    
    
    if (do_si_orth && is.null(oracle_si)) {
      
      joint_fit <- joint_orthogonalize_score(
        f_raw = f_use,
        X = X_use,
        u = u_use,
        theta_score = theta_score_use,
        fold_id = fold_use,
        N_si = test_fit$N_si,
        df_fallback = df_si_orth,
        C_lambda_beta = C_lambda_beta,
        C_lambda_theta = C_lambda_theta,
        eps = eps_num
      )
      
      f_orth_joint <- joint_fit$f_orth   # joint
      
      spline_orth_before <- joint_fit$spline_orth_before
      spline_orth_after <- joint_fit$spline_orth_after
      
      theta_orth_before <- joint_fit$theta_orth_before
      theta_orth_after <- joint_fit$theta_orth_after
      
      linear_orth_before_final <- joint_fit$linear_orth_before
      linear_orth_after_final <- joint_fit$linear_orth_after
      
      path_index <- joint_fit$path_index
      path_ratio <- joint_fit$path_ratio
      path_gap <- joint_fit$path_gap
      linear_kkt_ratio <- joint_fit$linear_kkt_ratio
      theta_kkt_ratio <- joint_fit$theta_kkt_ratio
      theta_kkt_ratio_by_fold <- joint_fit$theta_kkt_ratio_by_fold
      
      
    } else {
      
      # oracle / 关闭 SI 正交化时，仅线性正交化
      linear_final_fit <- linear_orthogonalize(f_use)
      f_orth_joint <- linear_final_fit$f_orth
      
      linear_orth_before_final <- linear_final_fit$linear_orth_before
      linear_orth_after_final <- linear_final_fit$linear_orth_after
    }
    
    
    ############################################################
    ## 7. 仅线性正交化结果，用于对比
    ############################################################
    
    linear_only_fit <- linear_orthogonalize(f_use)
    f_orth_linear <- linear_only_fit$f_orth
    
    
    ############################################################
    ## 8. 计算统计
    ############################################################
    
    calc_stat <- function(f_score) {
      
      rw_use <- r_use * f_score
      
      S_use <- sum(rw_use) / sqrt(n_use)
      V_use <- mean(rw_use^2)
      
      Zn_use <- S_use / sqrt(V_use + eps_num)
      p_use_value <- 1 - pnorm(Zn_use)
      
      list(
        Zn = Zn_use,
        p_value = p_use_value
      )
    }
    
    
    stat_joint <- calc_stat(f_orth_joint) 
    Zn <- stat_joint$Zn
    p_value <- stat_joint$p_value
    
    stat_linear <- calc_stat(f_orth_linear)
    Zn_linear <- stat_linear$Zn
    p_linear <- stat_linear$p_value
    
    
    ############################################################
    ## 9. 返回结果
    ############################################################
    
    list(
      
      Zn = Zn,
      p_value = p_value,
      n_use = n_use,
      
      Zn_linear = Zn_linear,
      p_linear = p_linear,
      
      # 线性正交化诊断
      linear_orth_before_final = linear_orth_before_final,
      linear_orth_after_final = linear_orth_after_final,
      
      # 单指标正交化诊断
      spline_orth_before = spline_orth_before,
      spline_orth_after = spline_orth_after,
      
      theta_orth_before = theta_orth_before,
      theta_orth_after = theta_orth_after,
      
      # KKT 
      path_index = path_index,
      path_ratio = path_ratio,
      path_gap = path_gap,
      linear_kkt_ratio = linear_kkt_ratio,
      theta_kkt_ratio = theta_kkt_ratio,
      theta_kkt_ratio_by_fold = theta_kkt_ratio_by_fold
    )
  }
  
  
  ############################################################
  ## 4. 两个方向
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
    if (any(!is.finite(pvals))) return(NA_real_)
    
    if (length(pvals) == 0) return(NA_real_)
    
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
  
  p_CCT_linear <- cct_pvalue(
    pvals = c(stat_12$p_linear, stat_21$p_linear)
  )
  
  
  ############################################################
  ## 6. 输出
  ############################################################
  
  list(
    p_CCT = p_CCT,
    p_CCT_linear = p_CCT_linear,
    
    id_part1 = id_part1,
    id_part2 = id_part2,
    
    fit_part1 = fit_part1,
    fit_part2 = fit_part2,
    
    stat_12 = stat_12,
    stat_21 = stat_21
  )
}



generate_PLSIM_data <- function(
    n = 1000, p = 600, q = 600, s = 20, c = 1,
    scenario = c("H0", "case1", "case2", "case3"),
    sigma0 = sqrt(1/2)
){
  scenario <- match.arg(scenario)
  
  # 1) 生成协变量
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  Z <- matrix(
    qnorm(runif(n * q, pnorm(-2.5), pnorm(2.5))),
    nrow = n,
    ncol = q
  )
  
  # 2) 单位向量系数
  beta <- numeric(p)               # 全部置 0
  beta[1:s] <- rnorm(s)            # 生成 s 个非零
  beta[1:s] <- beta[1:s] / sqrt(sum(beta[1:s]^2))
  
  theta <- numeric(q)                 # 全部置0
  theta[(q-s+1):q] <- rep(1,s)        # 后s个元素生成非零值
  theta <- theta / sqrt(sum(theta^2)) # 单位化
  
  # 3) 单指标 
  t_val  <- as.vector(Z %*% theta)
  linear_val = as.vector(X %*% beta)
  # Y_base <- as.vector(X %*% beta) + t_val^2  # data1
  # Y_base <- as.vector(X %*% beta) + cos(2*t_val)  # data2
  Y_base <- as.vector(X %*% beta) + exp(-t_val^2)  # data3
  
  dev <- rep(0, n)    # 默认无偏离，满足H0
  sigma <- rep(sigma0, n) 
  
  # 4) H1 场景
  if (scenario == "case1") {
    if (p < 6) stop("case1 需要 p >= 6")
    # t_val^2
    dev <- Z[,1]^3 + Z[,2]^2 + 0.5*exp(Z[,3])^3 + 3*abs(Z[,4]) + cos(Z[,5]+Z[,6])  
    
  } else if (scenario == "case2") {
    if (q < 6) stop("case2 需要 q >= 6")
    # cos(2*t_val)
    dev <- Z[,1]^3 + Z[,2]^2 + 2*exp(Z[,3])^3 + sin(Z[,4])^3 + 3*abs(Z[,5]) + abs(Z[,6])^3  
    
  } else if (scenario == "case3") {
    if (p < 7 || q < 3) stop("case3 需要 p >= 7 且 q >= 3")
    # exp(-t_val^2)
    dev <- 2*Z[,1]^3*abs(X[,1]) + Z[,2]^2 + exp(Z[,3])^3 + abs(Z[,3])*X[,2]^2
    
  } 
  
  # 5) 生成响应
  eps <- rnorm(n, mean = 0, sd = sigma)
  Y   <- Y_base + c*dev + eps
  
  list(Y = as.vector(Y), X = X, Z = Z, beta = beta, theta = theta, dev = dev, eps = eps)
}



#######################################################################################
start_time <- Sys.time()

n = 1000
p = 600
q = 600
s = 20

data <- generate_PLSIM_data(n = n, p = p, q = q, s = s, scenario = "H0", c = 1)
Y = data[[1]]
X = data[[2]]
Z = data[[3]]

result = PLSIM_test_split_ratio_cross2(
  Y=Y, X=X, Z=Z, 
  split_ratio = c(4, 4),
  K_si = 4,
  C_lambda_beta = 0.3,
  C_lambda_theta = 1,
  seed = 1
)

end_time = Sys.time()
end_time - start_time  # 15.87349 secs

# 最终两步正交化结果
cat("joint: \n")
result$stat_12$Zn
result$stat_12$p_value

result$stat_21$Zn
result$stat_21$p_value

result$p_CCT


cat("linear: \n")
result$stat_12$Zn_linear
result$stat_12$p_linear

result$stat_21$Zn_linear
result$stat_21$p_linear

result$p_CCT_linear




# =========================
# 并行计算：PLSIM 拟合优度检验的 size / power
# =========================
library(doParallel)
library(foreach)
library(doRNG)

# ====== Python 解释器与脚本路径 ======
py_bin  <- "D:/Anaconda_envs/envs/myenv/python.exe"  # change this path to your local Python environment.
py_file <- "code/python_code.py"

.EXPORT_FUNS <- c(
  "generate_PLSIM_data", 
  "build_si_basis", "build_oof_theta_score", 
  "joint_orthogonalize_score", "PLSIM_test_split_ratio_cross2",
  # 单指标相关函数
  ".safe_inv_sym", "Fd", "dot_Fd", "calc_N", "normalize_theta",
  "build_spline_matrices", "compute_dot_Bp", "compute_R_hat", "compute_hat_S_star",
  ".project_theta_minus_d", ".make_safe_gr", "estimate_theta", "estimate_g",
  "spline_single_index",  "predict_si", "set_R_python_seeds"
)


# -------------------------
# 统一设置 R + Python seed
# -------------------------
set_R_python_seeds <- function(seed, use_python = TRUE) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  
  if (use_python) {
    reticulate::py_run_string(sprintf("
import os, random
import numpy as np
try:
    import torch
except Exception:
    torch = None

seed = %d
os.environ['PYTHONHASHSEED'] = str(seed)
random.seed(seed)
np.random.seed(seed)

if torch is not None:
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    try:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    except Exception:
        pass
", as.integer(seed)))
  }
}


# -------------------------
# 并行功效 / 尺寸评估
# -------------------------
run_power_PLSIM <- function(
    n_sims   = 1000,
    n        = 1000,
    p        = 600,
    q        = 600,
    s        = 20,
    scenario = "H0",
    c        = 1,
    split_ratio = c(4, 4),
    K_si     = 4,
    C_lambda_beta = 0.3,
    C_lambda_theta = 1,
    ncores   = 10,
    seed     = 2026
) {
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  
  suppressPackageStartupMessages(library(reticulate))
  Sys.setenv(RETICULATE_PYTHON = py_bin)
  reticulate::use_python(py_bin, required = TRUE)
  reticulate::source_python(py_file)
  
  set_R_python_seeds(seed)
  
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  registerDoParallel(cl)
  
  parallel::clusterExport(
    cl,
    c("py_bin", "py_file", "seed", "set_R_python_seeds"),
    envir = environment()
  )
  
  parallel::clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1",
      NUMEXPR_NUM_THREADS = "1"
    )
    
    library(reticulate)
    Sys.setenv(RETICULATE_PYTHON = py_bin)
    reticulate::use_python(py_bin, required = TRUE)
    reticulate::source_python(py_file)
    
    set_R_python_seeds(seed)
    NULL
  })
  
  res_mat <- foreach(
    sim = 1:n_sims,
    .combine = rbind,
    .packages = c(
      "splines2", "MASS", "caret", "randomForest",
      "glmnet", "reticulate"
    ),
    .export = .EXPORT_FUNS,
    .options.RNG = seed
  ) %dorng% {
    
    set_R_python_seeds(seed + sim)
    
    tryCatch({
      
      dat <- generate_PLSIM_data(
        n = n, p = p, q = q, s = s, scenario = scenario, 
      )
      
      test_res <- PLSIM_test_split_ratio_cross2(
        Y = dat$Y, X = dat$X, Z = dat$Z, split_ratio = split_ratio, 
        K_si = K_si, C_lambda_beta = C_lambda_beta, C_lambda_theta = C_lambda_theta, 
        seed = seed + sim
      )
      
      c(
        p_12 = as.numeric(test_res$stat_12$p_value),
        p_21 = as.numeric(test_res$stat_21$p_value),
        p_CCT = as.numeric(test_res$p_CCT)
      )
      
    }, error = function(e) {
      
      warning(sprintf("sim %d failed: %s", sim, conditionMessage(e)))
      
      c(
        p_12 = NA_real_,
        p_21 = NA_real_,
        p_CCT = NA_real_
      )
    })
  }
  
  
  list(
    settings = list(
      n_sims = n_sims, n = n, p = p, q = q, s = s, scenario = scenario, 
      split_ratio = split_ratio, K_si = K_si, C_lambda_beta = C_lambda_beta,
      C_lambda_theta = C_lambda_theta, ncores = ncores, seed = seed
    ),
    
    raw_results = as.data.frame(res_mat)
    
  )
}


# =========================
# 使用示例
# =========================
start_time <- Sys.time()

res_size <- run_power_PLSIM(
  n_sims = 1000, n = 1000, p = 600, q = 600, s = 20, scenario = "H0", 
  split_ratio = c(4, 4), K_si = 4, C_lambda_beta = 0.3, C_lambda_theta = 1,
  ncores = 10, seed = 123
)


save(res_size, file = "results/high/exp(-t_val^2)_n1000_d600_s20_nsims1000_res_size.RData")

end_time <- Sys.time()
end_time - start_time  

rr = res_size$raw_results

# size/power
mean(rr$p_12 < 0.05, na.rm = TRUE)
mean(rr$p_21 < 0.05, na.rm = TRUE)
mean(rr$p_CCT < 0.05, na.rm = TRUE)

