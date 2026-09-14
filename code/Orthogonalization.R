# 构造样条基及其导数
build_si_basis <- function(u, N, derivs = 0L, eps = 1e-12) {
  u <- pmin(pmax(as.numeric(u), eps), 1 - eps)
  N <- as.integer(N)
  
  if (!is.finite(N) || N < 1L) {
    stop("N must be a positive integer.")
  }
  
  knots_si <- seq_len(N) / (N + 1)  
  
  if (derivs == 0L) {  # derivs=0 返回普通 B-spline 基
    out <- splines2::bSpline(
      x = u,
      knots = knots_si,
      degree = 3,
      Boundary.knots = c(0, 1),
      intercept = TRUE
    )
  } else {  # derivs=1 返回 B-spline 基关于u的一阶导数
    out <- splines2::dbs(
      x = u,
      derivs = as.integer(derivs),
      knots = knots_si,
      degree = 3,
      Boundary.knots = c(0, 1),
      intercept = TRUE
    )
  }
  
  as.matrix(out)
}


# 构造 θ 估计误差的一阶影响方向
build_oof_theta_score <- function(Z_new, fit_si, eps = 1e-12) {
  
  Z_new <- as.matrix(Z_new)
  
  required_names <- c(
    "theta_hat", "g_coef_hat", "Z_mu", "Z_sd", "N", "a")
  
  missing_names <- setdiff(required_names, names(fit_si))
  if (length(missing_names) > 0L) {
    stop(
      "fit_si is missing: ",
      paste(missing_names, collapse = ", ")
    )
  }
  
  Z_sd_use <- as.numeric(fit_si$Z_sd)
  Z_sd_use[!is.finite(Z_sd_use) | Z_sd_use == 0] <- 1
  
  Z_std <- sweep(
    sweep(Z_new, MARGIN = 2, STATS = as.numeric(fit_si$Z_mu), FUN = "-"),
    MARGIN = 2, STATS = Z_sd_use, FUN = "/"
  )
  
  theta_hat <- as.numeric(fit_si$theta_hat)
  theta_norm <- sqrt(sum(theta_hat^2))
  
  if (!is.finite(theta_norm) || theta_norm <= eps) {
    stop("Invalid theta_hat.")
  }
  
  theta_hat <- theta_hat / theta_norm
  
  if (length(theta_hat) != ncol(Z_std)) {
    stop("The dimensions of theta_hat and Z_new do not match.")
  }
  
  v <- as.vector(Z_std %*% theta_hat)
  u <- Fd(v, d = ncol(Z_std), a = fit_si$a, eps = eps)
  
  dB_du <- build_si_basis(
    u = u,
    N = fit_si$N,
    derivs = 1L,
    eps = eps
  )
  
  g_coef <- as.numeric(fit_si$g_coef_hat)
  
  if (ncol(dB_du) != length(g_coef)) {
    stop(
      "The spline derivative basis has ",
      ncol(dB_du),
      " columns, but g_coef_hat has length ",
      length(g_coef),
      "."
    )
  }
  
  g_prime_u <- as.vector(dB_du %*% g_coef)
  F_prime_v <- dot_Fd(v, d = ncol(Z_std), a = fit_si$a)
  
  # theta 的单位球切空间
  P_theta <- diag(length(theta_hat)) - tcrossprod(theta_hat)
  theta_tangent <- Z_std %*% P_theta
  
  theta_score <- sweep(theta_tangent, MARGIN = 1, STATS = g_prime_u * F_prime_v, FUN = "*")
  
  theta_score[!is.finite(theta_score)] <- 0
  
  list(
    u = pmin(pmax(u, eps), 1 - eps),
    theta_score = theta_score
  )
}


# ============================================================
# partially penalized square-root Lasso
# ============================================================
joint_orthogonalize_score <- function(
    f_raw, X, u, theta_score, fold_id, N_si = NULL, df_fallback = 6,
    C_lambda_beta = 0.3, C_lambda_theta = 1, eps = 1e-10) 
{
  
  f_raw <- as.numeric(f_raw)
  X <- as.matrix(X)
  u <- as.numeric(u)
  theta_score <- as.matrix(theta_score)
  fold_id <- as.integer(fold_id)
  
  n <- length(f_raw)
  p <- ncol(X)
  q <- ncol(theta_score)
  
  stopifnot(nrow(X) == n, length(u) == n, nrow(theta_score) == n, length(fold_id) == n)
  
  if (any(!is.finite(f_raw)) || any(!is.finite(X)) ||
      any(!is.finite(u)) || any(!is.finite(theta_score))) {
    stop("Non-finite values are not allowed in joint orthogonalization.")
  }
  
  fold_levels <- sort(unique(fold_id))
  K_use <- length(fold_levels)
  
  N_use <- if (!is.null(N_si) && is.finite(N_si)) {
    max(as.integer(N_si), 1L)
  } else {
    max(as.integer(df_fallback), 1L)
  }
  
  ############################################################
  ## 1. fold-specific spline blocks
  ############################################################
  
  B_design <- matrix(numeric(0), nrow = n, ncol = 0)
  B_info <- vector("list", K_use)
  
  for (kk in seq_along(fold_levels)) {
    
    k <- fold_levels[kk]
    idx <- which(fold_id == k)
    
    Bk_full <- build_si_basis(
      u = u[idx],
      N = N_use,
      derivs = 0L,
      eps = eps
    )
    
    qr_B <- qr(Bk_full, tol = 1e-8)
    if (qr_B$rank > 0L) {
      keep_B <- qr_B$pivot[seq_len(qr_B$rank)]
      Bk <- Bk_full[, keep_B, drop = FALSE]
    } else {
      keep_B <- integer(0)
      Bk <- matrix(numeric(0), nrow = length(idx), ncol = 0)
    }
    
    start_col <- ncol(B_design) + 1L
    
    if (ncol(Bk) > 0L) {
      block <- matrix(0, nrow = n, ncol = ncol(Bk))
      block[idx, ] <- Bk
      B_design <- cbind(B_design, block)
      cols <- start_col:(start_col + ncol(Bk) - 1L)
    } else {
      cols <- integer(0)
    }
    
    B_info[[kk]] <- list(
      fold = k,
      idx = idx,
      B = Bk,
      keep = keep_B,
      cols = cols
    )
  }
  
  
  ############################################################
  ## 2. fold-specific theta-score blocks
  ############################################################
  
  S_design <- matrix(numeric(0), nrow = n, ncol = 0)
  S_info <- vector("list", K_use)
  theta_penalty_raw <- numeric(0)
  lambda_theta_vec <- numeric(K_use)
  
  for (kk in seq_along(fold_levels)) {
    
    k <- fold_levels[kk]
    idx <- which(fold_id == k)
    nk <- length(idx)
    
    Sk_full <- theta_score[idx, , drop = FALSE]
    
    finite_col <- apply(Sk_full, 2, function(v) all(is.finite(v)))
    scale_col <- apply(Sk_full, 2, function(v) sqrt(mean(v^2)))
    keep_S <- finite_col & is.finite(scale_col) & scale_col > eps
    
    Sk <- Sk_full[, keep_S, drop = FALSE]
    
    theta_scale <- scale_col[keep_S]
    
    lambda_theta_k <- C_lambda_theta * sqrt(log(max(q, 2L)) / nk)
    lambda_theta_vec[kk] <- lambda_theta_k
    
    start_col <- ncol(S_design) + 1L
    
    if (ncol(Sk) > 0L) {
      block <- matrix(0, nrow = n, ncol = ncol(Sk))
      block[idx, ] <- Sk
      S_design <- cbind(S_design, block)
      cols <- start_col:(start_col + ncol(Sk) - 1L)
      
      theta_penalty_raw <- c(theta_penalty_raw, (nk / n) * lambda_theta_k * theta_scale)
      
    } else {
      cols <- integer(0)
    }
    
    S_info[[kk]] <- list(
      fold = k,
      idx = idx,
      S = Sk,
      keep = which(keep_S),
      cols = cols,
      lambda = lambda_theta_k,
      scale = theta_scale       
    )
  }
  
  
  ############################################################
  ## 3. 构造联合 design： [B_1,...,B_K, X, S_1,...,S_K]
  ############################################################
  
  D_joint <- cbind(B_design, X, S_design)
  
  n_B <- ncol(B_design)
  n_X <- ncol(X)
  n_S <- ncol(S_design)
  
  lambda_beta <- C_lambda_beta * sqrt(log(max(p, 2L)) / n)
  
  penalty_raw <- c(rep(0, n_B), rep(lambda_beta, n_X), theta_penalty_raw)
  
  
  ############################################################
  ## 4. Square-root Lasso via the complete Lasso path
  ############################################################
  
  if (sum(penalty_raw) > 0) {
    
    # glmnet 会把 penalty.factor 重标到和为变量总数
    pf_scale <- ncol(D_joint) / sum(penalty_raw)
    penalty_factor <- pf_scale * penalty_raw
    
    fit_joint <- glmnet::glmnet(
      x = D_joint,
      y = f_raw,
      alpha = 1,
      intercept = FALSE,
      standardize = FALSE,
      penalty.factor = penalty_factor,
      nlambda = 1000,
      thresh = 1e-10,
      maxit = 1000000
    )
    
    pred_path <- as.matrix(D_joint %*% fit_joint$beta)
    
    resid_path <- sweep(pred_path, 1, f_raw, function(pred, obs) obs - pred)
    
    scale_path <- sqrt(colMeans(resid_path^2))
    scale_path <- pmax(scale_path, eps)
    
    # Square-root Lasso matching condition: lambda_glmnet * pf_scale ≈ residual scale
    index_joint <- which.min(
      abs(fit_joint$lambda * pf_scale - scale_path)
    )
    
    coef_joint <- as.numeric(fit_joint$beta[, index_joint])
    
    f_orth <- as.numeric(resid_path[, index_joint])
    
    sigma_hat <- scale_path[index_joint]
    lambda_glmnet_selected <- fit_joint$lambda[index_joint]
    lambda_target_selected <- sigma_hat / pf_scale
    
    path_ratio <- (lambda_glmnet_selected * pf_scale) / sigma_hat  
    path_gap <- abs(lambda_glmnet_selected * pf_scale - sigma_hat) 
    
  } else {
    
    fit_ls <- lm.fit(
      x = D_joint,
      y = f_raw
    )
    
    coef_joint <- as.numeric(fit_ls$coefficients)
    coef_joint[!is.finite(coef_joint)] <- 0
    
    f_orth <- f_raw - as.vector(D_joint %*% coef_joint)
    
    sigma_hat <- sqrt(mean(f_orth^2))
    
    pf_scale <- NA_real_
    index_joint <- NA_integer_
    lambda_glmnet_selected <- NA_real_
    lambda_target_selected <- NA_real_
    path_ratio <- NA_real_
    path_gap <- NA_real_
  }
  
  
  ############################################################
  ## 5. 提取各系数块
  ############################################################
  
  eta_hat <- vector("list", K_use)
  for (kk in seq_along(B_info)) {
    cols <- B_info[[kk]]$cols
    eta_hat[[kk]] <- if (length(cols) > 0L) {
      coef_joint[cols]
    } else {
      numeric(0)
    }
  }
  
  
  beta_cols <- if (n_X > 0L) {
    (n_B + 1L):(n_B + n_X)
  } else {
    integer(0)
  }
  
  omega_beta_hat <- if (length(beta_cols) > 0L) {
    coef_joint[beta_cols]
  } else {
    numeric(0)
  }
  
  theta_eta_hat <- vector("list", K_use)
  theta_selected <- 0L
  
  if (n_S > 0L) {
    S_offset <- n_B + n_X
    
    for (kk in seq_along(S_info)) {
      cols_local <- S_info[[kk]]$cols
      
      coef_full <- rep(0, q)
      
      if (length(cols_local) > 0L) {
        coef_kept <- coef_joint[S_offset + cols_local]
        coef_full[S_info[[kk]]$keep] <- coef_kept
        theta_selected <- theta_selected + sum(abs(coef_kept) > 1e-8)
      }
      
      theta_eta_hat[[kk]] <- coef_full
    }
  } else {
    theta_eta_hat <- lapply(seq_len(K_use), function(.) rep(0, q))
  }
  
  
  ############################################################
  ## 6. KKT / orthogonality diagnostics
  ############################################################
  
  spline_before_vec <- spline_after_vec <- numeric(K_use)
  theta_before_vec <- theta_after_vec <- numeric(K_use)
  
  for (kk in seq_along(fold_levels)) {
    
    idx <- B_info[[kk]]$idx
    nk <- length(idx)
    
    Bk <- B_info[[kk]]$B
    spline_before_vec[kk] <- if (ncol(Bk) > 0L) {
      max(abs(as.vector(crossprod(Bk, f_raw[idx])) / nk))
    } else {
      0
    }
    spline_after_vec[kk] <- if (ncol(Bk) > 0L) {
      max(abs(as.vector(crossprod(Bk, f_orth[idx])) / nk))
    } else {
      0
    }
    
    Sk <- S_info[[kk]]$S
    theta_before_vec[kk] <- if (ncol(Sk) > 0L) {
      max(abs(as.vector(crossprod(Sk, f_raw[idx])) / nk))
    } else {
      0
    }
    theta_after_vec[kk] <- if (ncol(Sk) > 0L) {
      max(abs(as.vector(crossprod(Sk, f_orth[idx])) / nk))
    } else {
      0
    }
  }
  
  linear_orth_before <- if (n_X > 0L) {
    max(abs(as.vector(crossprod(X, f_raw)) / n))
  } else {
    0
  }
  
  linear_orth_after <- if (n_X > 0L) {
    max(abs(as.vector(crossprod(X, f_orth)) / n))
  } else {
    0
  }
  
  # KKT ratio diagnostics
  penalty_scale_selected <- lambda_glmnet_selected * pf_scale
  
  # X block：保持原来的 KKT
  linear_kkt_ratio <- linear_orth_after / (penalty_scale_selected * lambda_beta)
  
  # theta block：考虑每一列 theta-score 的 scale loading
  theta_kkt_ratio_by_fold <- numeric(K_use)
  
  for (kk in seq_along(fold_levels)) {
    
    idx <- S_info[[kk]]$idx
    nk <- length(idx)
    Sk <- S_info[[kk]]$S
    theta_scale <- S_info[[kk]]$scale
    
    if (ncol(Sk) > 0L) {
      
      theta_grad <- abs(
        as.vector(crossprod(Sk, f_orth[idx])) / nk
      )
      
      theta_kkt_ratio_by_fold[kk] <- max(
        theta_grad / (penalty_scale_selected * lambda_theta_vec[kk] * theta_scale)
      )
      
    } else {
      
      theta_kkt_ratio_by_fold[kk] <- 0
      
    }
  }
  
  theta_kkt_ratio <- max(theta_kkt_ratio_by_fold)
  
  
  list(
    f_orth = f_orth,
    
    spline_orth_before = max(spline_before_vec),
    spline_orth_after = max(spline_after_vec),
    
    linear_orth_before = linear_orth_before,
    linear_orth_after = linear_orth_after,
    
    theta_orth_before = max(theta_before_vec),
    theta_orth_after = max(theta_after_vec),

    path_index = index_joint,
    path_ratio = path_ratio,
    path_gap = path_gap,
    linear_kkt_ratio = linear_kkt_ratio,
    theta_kkt_ratio = theta_kkt_ratio,
    theta_kkt_ratio_by_fold = theta_kkt_ratio_by_fold
  )
}

