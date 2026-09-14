# ============================================================
# 基于之前使用的 SIM_functions.R 修改，解决高维速度过慢的问题。
# 主要修改：
# 1) compute_R_hat 不再显式构造 n×n 投影矩阵；
# 2) compute_hat_S_star 不再显式构造 P_theta，且 dotFd 使用稳定版 dot_Fd；
# 3) estimate_theta 在高维时仍用同一准则/同一梯度，但减少多启动和兜底优化。
# ============================================================

library(splines2)
library(MASS)
library(caret)

# ========================
# 稳健地求一个对称矩阵的逆
# ========================
.safe_inv_sym <- function(M, ridge = 0, max_boost = 5L) {
  p <- ncol(M)
  lam <- ridge
  
  for (k in 0:max_boost) {
    Mi <- M + (lam * diag(p))
    ch <- try(chol(Mi), silent = TRUE)
    if (!inherits(ch, "try-error")) {
      return(list(inv = chol2inv(ch), ridge_used = lam, method = "chol"))
    }
    lam <- if (lam == 0) 1e-8 else lam * 10
  }
  
  sv <- svd(M)
  tol <- max(dim(M)) * .Machine$double.eps * max(sv$d, 1)
  d_inv <- ifelse(sv$d > tol, 1 / sv$d, 0)
  inv <- sv$u %*% (d_inv * t(sv$v))
  list(inv = inv, ridge_used = lam, method = "svd")
}


# ========================
# Fd 与其导数
# ========================
Fd <- function(nu, d, a, eps = 1e-8) {
  nu <- pmin(pmax(nu, -a), a)
  shape <- (d + 1) / 2
  t <- nu / a
  U <- pbeta((t + 1) / 2, shape1 = shape, shape2 = shape)
  U <- pmin(pmax(U, eps), 1 - eps)
  return(U)
}

## 容易得到NaN
# dot_Fd <- function(x, d, a) {
#   if (abs(x) > a) return(0)
#   shape <- (d + 1) / 2
#   numerator <- gamma(d + 1)
#   denominator <- a * (gamma(shape))^2 * (2^d)
#   term <- (1 - (x / a)^2)^((d - 1) / 2)
#   numerator / denominator * term
# }


dot_Fd <- function(x, d, a) {
  shape <- (d + 1) / 2
  out <- numeric(length(x))
  mask <- abs(x) <= a
  
  z <- x[mask] / a
  one_minus_z2 <- pmax(1 - z^2, 0)
  
  # 大 d 时用对数版本，因为gamma(d + 1)和gamma(shape)^2会溢出，变成InF
  if (d > 50) {
    log_const <- lgamma(d + 1) - lgamma(shape) - lgamma(shape) - d * log(2) - log(a)
    log_val <- log_const + ((d - 1) / 2) * log(one_minus_z2)
    out[mask] <- exp(log_val)
  } else {
    # 小 d 时用原版（更快）
    out[mask] <- gamma(d + 1) / (a * gamma(shape)^2 * 2^d) * (one_minus_z2)^((d - 1) / 2)
  }
  
  out
}


# ========================
# N 计算
# ========================
calc_N <- function(n, c1 = 1, c2 = 5) {
  N1 <- c1 * floor(n^(1 / 5.5))
  N <- min(N1, c2)
  return(max(N, 4))  # 之前是max(N, 2)
}


# ========================
# theta 规范化
# ========================
normalize_theta <- function(theta) {
  theta_norm <- theta / sqrt(sum(theta^2))
  d <- length(theta_norm)
  if (theta_norm[d] < 0) theta_norm <- -theta_norm
  theta_norm
}


# ========================
#  B-spline 样条矩阵与投影矩阵
# ========================
build_spline_matrices <- function(U_theta, N) {
  knots <- seq(from = 1/(N+1), to = N/(N+1), by = 1/(N+1))
  B_theta <- bSpline(
    x = U_theta,
    knots = knots,
    degree = 3,
    Boundary.knots = c(0, 1),
    intercept = TRUE
  )
  BtB <- t(B_theta) %*% B_theta
  inv_res <- .safe_inv_sym(BtB, ridge = getOption("SI_ridge", 0))
  P_theta <- B_theta %*% inv_res$inv %*% t(B_theta)
  list(B_theta = as.matrix(B_theta), P_theta = P_theta, BtB_inv = inv_res$inv)
}


# ========================
# 计算 dB/dθ_p
# ========================
compute_dot_Bp <- function(U_theta, X_std, theta, p, N, a) {
  n <- nrow(X_std)
  d <- ncol(X_std)
  
  X_theta  <- as.vector(X_std %*% theta)
  fd_prime <- vapply(X_theta, function(x) dot_Fd(x, d = d, a = a), numeric(1))
  Xp       <- X_std[, p]
  
  knots <- seq(1 / (N + 1), N / (N + 1), by = 1 / (N + 1))
  dB3_dU <- as.matrix(splines2::dbs(
    x = U_theta,
    knots = knots,
    degree = 3,
    Boundary.knots = c(0, 1),
    derivs = 1,
    intercept = TRUE
  ))
  dot_Bp <- dB3_dU * (fd_prime %o% rep(1, ncol(dB3_dU)))
  dot_Bp <- dot_Bp * (Xp %o% rep(1, ncol(dB3_dU)))
  dot_Bp
}


# ========================
# 经验风险（目标函数）：R_hat^*
# 改动：不再显式构造 n×n 的 P_theta，而是链式计算 P_theta Y。
# 数学等价于 B(B'B)^(-1)B'Y，但高维/大样本时快很多。
# ========================
compute_R_hat <- function(theta_minus_d, X_std, Y, a, N, ridge = getOption("SI_ridge", 0)) {
  d <- ncol(X_std)
  nsq <- sum(theta_minus_d^2)
  if (nsq >= 1) theta_minus_d <- theta_minus_d / sqrt(nsq) * 0.999
  
  theta_d <- sqrt(max(0, 1 - sum(theta_minus_d^2)))
  theta <- normalize_theta(c(theta_minus_d, theta_d))
  
  X_theta <- as.vector(X_std %*% theta)
  U_theta <- pmin(pmax(Fd(X_theta, d = d, a = a), 1e-12), 1 - 1e-12)
  
  knots <- seq(1/(N+1), N/(N+1), by = 1/(N+1))
  B_theta <- as.matrix(splines2::bSpline(
    x = U_theta,
    knots = knots,
    degree = 3,
    Boundary.knots = c(0, 1),
    intercept = TRUE
  ))
  
  BtB <- crossprod(B_theta)
  inv_res <- .safe_inv_sym(BtB, ridge = ridge)
  Y_hat <- as.vector(B_theta %*% inv_res$inv %*% crossprod(B_theta, Y))
  
  mean((Y - Y_hat)^2)  # profile least squares
}


# ========================
# score vector S_hat^*，对应论文 Lemma 3.1
# 改动：
# 1) 不显式构造 n×n 的 P_theta；
# 2) dotFd 统一调用稳定版 dot_Fd()，高维时用 lgamma 避免溢出。
# ========================
compute_hat_S_star <- function(theta_minus_d, X_std, Y, a, N, ridge = getOption("SI_ridge", 0)) {
  d <- ncol(X_std); n <- nrow(X_std)
  
  # 约束：||theta_-d|| < 1，避免 theta_d→0
  nsq <- sum(theta_minus_d^2)
  if (nsq >= 1) theta_minus_d <- theta_minus_d / sqrt(nsq) * 0.999
  
  # 还原并标准化到上半球
  theta_d <- sqrt(max(0, 1 - sum(theta_minus_d^2)))
  theta <- normalize_theta(c(theta_minus_d, theta_d))
  
  # 单指标与U
  X_theta <- as.vector(X_std %*% theta)
  U_theta <- pmin(pmax(Fd(X_theta, d = d, a = a), 1e-12), 1 - 1e-12)
  
  # 样条矩阵/逆：一次构建，后面共用
  knots <- seq(1/(N+1), N/(N+1), by = 1/(N+1))
  B_theta <- as.matrix(splines2::bSpline(
    x = U_theta, knots = knots, degree = 3,
    Boundary.knots = c(0, 1), intercept = TRUE
  ))
  BtB <- crossprod(B_theta)
  inv_res <- .safe_inv_sym(BtB, ridge = ridge)
  BtB_inv <- inv_res$inv
  
  # P_theta Y —— 链式计算，避免显式构造 P_theta
  BY <- crossprod(B_theta, Y)
  t2 <- BtB_inv %*% BY
  P_Y <- as.vector(B_theta %*% t2)
  Y_res <- as.vector(Y - P_Y)
  
  # dB/dU（导数基）
  dB_dU <- as.matrix(splines2::dbs(
    x = U_theta, knots = knots, degree = 3,
    Boundary.knots = c(0, 1), derivs = 1, intercept = TRUE
  ))
  
  # 高维稳定版 F'_d(Xθ)：不能直接用 gamma(d+1)
  dotFd <- dot_Fd(X_theta, d = d, a = a)
  
  # t3 = (dB/dU) (B'B)^(-1) B'Y
  t3 <- as.vector(dB_dU %*% t2)
  q <- Y_res * t3
  
  # d 分量
  Y_dot_Pd_Y_const <- 2 * sum(q * dotFd * X_std[, d])
  
  # 1,...,d-1 分量，一次性向量化
  Xp_mat <- X_std[, 1:(d-1), drop = FALSE]
  Y_dot_Pp_Y_vec <- 2 * colSums(Xp_mat * (dotFd * q))
  
  hat_S_star <- -(1 / n) * (
    Y_dot_Pp_Y_vec -
      (theta[1:(d-1)] / max(theta[d], 1e-12)) * Y_dot_Pd_Y_const
  )
  
  if (!all(is.finite(hat_S_star))) hat_S_star[] <- NA_real_
  hat_S_star
}


# ========================
# θ_-d 投影
# ========================
.project_theta_minus_d <- function(theta_minus_d, c_margin = 0.05) {
  r_max <- sqrt(max(1e-12, 1 - c_margin^2))
  nrm <- sqrt(sum(theta_minus_d^2))
  if (nrm >= r_max) theta_minus_d <- theta_minus_d / nrm * (r_max * 0.999)
  theta_minus_d
}


# ========================
# 安全梯度壳：解析优先，失败→有限差分
# ========================
.make_safe_gr <- function(objective_fun, analytic_gr, eps = 1e-6) {
  function(theta) {
    g <- try(analytic_gr(theta), silent = TRUE)
    if (!inherits(g, "try-error") && all(is.finite(g))) return(g)
    # 兜底：有限差分
    f0 <- objective_fun(theta); p <- length(theta); gd <- numeric(p)
    for (j in 1:p) {
      e <- rep(0, p); e[j] <- eps
      f1 <- objective_fun(theta + e); if (!is.finite(f1)) f1 <- f0
      gd[j] <- (f1 - f0) / eps
    }
    gd
  }
}


# ========================
# 估计 theta：复现论文 Step 3
# ========================
estimate_theta <- function(X_std, Y, a_si, 
                           c_margin = 0.01,
                           a_q = 0.995,
                           maxN = 5,
                           ridge = 1e-6,
                           n_starts = 5,
                           grad_tol = 1e-4,
                           rel_impr_tol = 1e-8,
                           # 高维加速控制
                           high_dim_fast = TRUE,
                           high_dim_threshold = 100,
                           high_dim_a_q = 0.995,
                           high_dim_n_starts = 1,
                           high_dim_eval_max = 500,
                           high_dim_iter_max = 500,
                           high_dim_use_fallback = TRUE,
                           high_dim_use_lse_init = FALSE,
                           verbose = FALSE) {
  
  X_std <- as.matrix(X_std)
  Y <- as.numeric(Y)
  d <- ncol(X_std)
  n <- nrow(X_std)
  
  is_high_dim <- high_dim_fast && (d > high_dim_threshold)
  
  # 高维仍然用同一套准则，只调整计算和优化控制。
  a_q_use <- if (is_high_dim) high_dim_a_q else a_q
  n_starts_use <- if (is_high_dim) high_dim_n_starts else n_starts
  eval_max_use <- if (is_high_dim) high_dim_eval_max else 1000
  iter_max_use <- if (is_high_dim) high_dim_iter_max else 1000
  use_fallback <- if (is_high_dim) high_dim_use_fallback else TRUE
  use_lse_init <- if (is_high_dim) high_dim_use_lse_init else TRUE
  
  if (verbose) {
    cat(sprintf(
      "SI theta: d=%d, n=%d, mode=%s, n_starts=%d, eval.max=%d, iter.max=%d\n",
      d, n, ifelse(is_high_dim, "high-dim accelerated", "low-dim robust"),
      n_starts_use, eval_max_use, iter_max_use
    ))
  }
  
  if (is.null(a_si)) {
    X_norm <- sqrt(rowSums(X_std^2))
    a <- as.numeric(quantile(X_norm, a_q_use))
  } else {
    a <- as.numeric(a_si)
  }
  
  N <- min(calc_N(n = n), maxN)
  
  .proj <- function(v) {
    r_max <- sqrt(max(1e-12, 1 - c_margin^2))
    nv <- sqrt(sum(v^2))
    if (nv >= r_max) v <- v / nv * (r_max * 0.999)
    v
  }
  
  objective_fun_raw <- function(theta_minus_d) {
    theta_minus_d <- .proj(theta_minus_d)
    old <- options(SI_ridge = ridge); on.exit(options(old), add = TRUE)
    compute_R_hat(theta_minus_d, X_std, Y, a, N, ridge = ridge)
  }
  analytic_gr_raw <- function(theta_minus_d) {
    theta_minus_d <- .proj(theta_minus_d)
    old <- options(SI_ridge = ridge); on.exit(options(old), add = TRUE)
    compute_hat_S_star(theta_minus_d, X_std, Y, a, N, ridge = ridge)
  }
  gradient_fun <- .make_safe_gr(objective_fun_raw, analytic_gr_raw)
  
  # 初值：低维保留 LSE 初值；高维默认不用 LSE，避免 lm 在 p>=n 或大 p 情形极慢/不稳。
  inits <- list()
  if (use_lse_init && d < n) {
    lse <- try(coef(lm(Y ~ X_std - 1)), silent = TRUE)
    if (!inherits(lse, "try-error") && all(is.finite(lse))) {
      lse <- as.numeric(lse)
      if (lse[d] < 0) lse <- -lse
      lse <- lse / sqrt(sum(lse^2))
      inits[[length(inits)+1]] <- .proj(lse[1:(d-1)])
    }
  }
  if (length(inits) == 0) {
    inits[[length(inits)+1]] <- .proj(rep(0, d - 1))
  }
  
  if (n_starts_use >= 2) {
    for (t in 2:n_starts_use) {
      v <- rnorm(d - 1)
      v <- v / sqrt(sum(v^2) + 1e-12) *
        runif(1, 0, sqrt(1 - c_margin^2) * 0.8)
      inits[[length(inits)+1]] <- .proj(v)
    }
  }
  
  best <- list(value = Inf, par = NULL, g_inf = Inf, conv_ok = FALSE,
               convergence = NA_integer_, message = NA_character_)
  
  one_try <- function(init_v) {
    fit <- nlminb(
      start = init_v,
      objective = objective_fun_raw,
      gradient  = gradient_fun,
      lower = rep(-Inf, d-1),
      upper = rep( Inf, d-1),
      control = list(
        eval.max = eval_max_use,
        iter.max = iter_max_use,
        rel.tol = 1e-8,
        x.tol = 1e-8
      )
    )
    
    par <- .proj(fit$par)
    val <- objective_fun_raw(par)
    g   <- gradient_fun(par)
    g_inf <- max(abs(g))
    conv_ok <- is.finite(g_inf) && (g_inf <= grad_tol || fit$convergence == 0)
    
    # nlminb() 做优化如果梯度不够小或者收敛不好，再用 optim(..., method="L-BFGS-B") 兜底
    if (use_fallback && !conv_ok) {
      fit2 <- optim(
        par    = par,
        fn     = objective_fun_raw,
        gr     = gradient_fun,
        method = "L-BFGS-B",
        lower  = rep(-0.999, d-1),
        upper  = rep( 0.999, d-1),
        control= list(maxit = eval_max_use, pgtol = 1e-8)
      )
      par <- .proj(fit2$par)
      val <- objective_fun_raw(par)
      g   <- gradient_fun(par)
      g_inf <- max(abs(g))
      conv_ok <- is.finite(g_inf) && (g_inf <= grad_tol || fit2$convergence == 0)
      return(list(par=par, value=val, g_inf=g_inf, conv_ok=conv_ok,
                  convergence=fit2$convergence, message="fallback optim"))
    }
    
    list(par=par, value=val, g_inf=g_inf, conv_ok=conv_ok,
         convergence=fit$convergence, message=fit$message)
  }
  
  for (z in inits) {
    res <- one_try(z)
    if (is.finite(res$value) && res$value < best$value) best <- res
  }
  
  if (is.null(best$par)) {
    # 极端情况下的兜底，避免函数直接报错。
    best$par <- .proj(rep(0, d - 1))
    best$value <- objective_fun_raw(best$par)
    best$g_inf <- max(abs(gradient_fun(best$par)))
    best$conv_ok <- FALSE
    best$convergence <- 99L
    best$message <- "all starts failed; returned zero start"
  }
  
  theta_minus_d_hat <- .proj(best$par)
  theta_d_hat <- sqrt(max(0, 1 - sum(theta_minus_d_hat^2)))
  theta_hat <- normalize_theta(c(theta_minus_d_hat, theta_d_hat))
  
  list(
    theta_hat = theta_hat,
    a = a,
    N = N,
    version = ifelse(is_high_dim, "same-criterion accelerated high-dim", "robust low-dim"),
    opt_message = if (best$conv_ok) "converged/accepted" else "accepted with large gradient",
    grad_inf = best$g_inf,
    value = best$value,
    convergence = best$convergence,
    message = best$message,
    is_high_dim = is_high_dim
  )
}





# ========================
# 估计 g：对应论文 Step 4 和公式 (2.10)
# ========================
estimate_g <- function(theta_hat, X_std, Y, a, N) {
  
  X_theta_hat <- as.vector(X_std %*% theta_hat)
  U_theta_hat <- Fd(X_theta_hat, d = ncol(X_std), a = a)
  U_theta_hat <- pmin(pmax(U_theta_hat, 1e-12), 1 - 1e-12)
  
  knots <- seq(1/(N+1), N/(N+1), by = 1/(N+1))
  B_theta_hat <- bSpline(
    x = U_theta_hat,
    knots = knots,
    degree = 3,
    Boundary.knots = c(0, 1),
    intercept = TRUE
  )
  BtB <- t(B_theta_hat) %*% B_theta_hat
  inv_res <- .safe_inv_sym(BtB, ridge = getOption("SI_ridge", 0))
  beta_hat <- inv_res$inv %*% t(B_theta_hat) %*% Y
  
  coef_hat <- inv_res$inv %*% crossprod(B_theta_hat, Y)  # spline least squares estimator
  
  g_hat <- function(v) {
    U <- pmin(pmax(Fd(v, d = ncol(X_std), a = a), 1e-12), 1 - 1e-12)
    B_U <- bSpline(
      x = U, knots = knots, degree = 3,
      Boundary.knots = c(0, 1), intercept = TRUE
    )
    as.vector(B_U %*% beta_hat)
  }
  
  # 数值导数 g'(v) —— 中心差分，v 可以是向量
  g_prime_hat <- function(v, h = 1e-4) {
    v <- as.numeric(v)
    (g_hat(v + h) - g_hat(v - h)) / (2 * h)
  }
  
  list(g_hat = g_hat, 
       g_prime_hat = g_prime_hat, 
       X_theta_hat = X_theta_hat,
       coef_hat = coef_hat   
       )
}



# ========================
# 主函数：输入原始 Z，内部标准化
# ========================
spline_single_index <- function(Z_train, Y, Z_mu, Z_sd, a_si) {
  
  Z_train <- as.matrix(Z_train)
  Y <- as.numeric(Y)
  Z_mu <- Z_mu 
  Z_sd <- Z_sd
  Z_std <- sweep(sweep(Z_train, 2, Z_mu, "-"), 2, Z_sd, "/")
  
  theta_result <- estimate_theta(Z_std, Y, a_si)
  
  g_result <- estimate_g(
    theta_hat = theta_result$theta_hat,
    X_std = Z_std,
    Y = Y,
    a = a_si,
    N = theta_result$N
  )
  
  Y_hat <- g_result$g_hat(g_result$X_theta_hat)
  mse <- mean((Y - Y_hat)^2)
  
  theta_hat_std <- theta_result$theta_hat  # 标准化尺度theta

  
  list(
    theta_hat = theta_hat_std,
    g_hat = g_result$g_hat,
    g_prime_hat = g_result$g_prime_hat, 
    X_theta_hat = g_result$X_theta_hat,
    a = theta_result$a,
    N = theta_result$N,
    Y_hat = Y_hat,
    mse = mse,
    
    g_coef_hat = g_result$coef_hat,   
    Z_mu = Z_mu,
    Z_sd = Z_sd,
    theta_version = theta_result$version,
    theta_value = theta_result$value,
    theta_grad_inf = theta_result$grad_inf,
    theta_convergence = theta_result$convergence,
    theta_message = theta_result$message,
    theta_is_high_dim = theta_result$is_high_dim
  )
}


## =====================================================
## predict：输入原始 Z_new，返回 R 的预测值
## =====================================================

predict_si <- function(theta, g_coef, Z_new, Z_mu, Z_sd, N, a, d) {
  
  Z_new <- as.matrix(Z_new)
  
  # 标准化（必须用训练时的 mu 和 sd）
  Z_new_std <- sweep(sweep(Z_new, 2, Z_mu, "-"), 2, Z_sd, "/")
  
  # 指标
  index_new <- as.vector(Z_new_std %*% theta)
  
  # 构造样条基（和 estimate_g 完全一样）
  knots <- seq(1/(N+1), N/(N+1), by = 1/(N+1))
  U <- Fd(index_new, d = d, a = a)
  B_new <- as.matrix(splines2::bSpline(
    x = U,
    knots = knots,
    degree = 3,
    Boundary.knots = c(0,1),
    intercept = TRUE
  ))
  
  as.vector(B_new %*% g_coef)
}


