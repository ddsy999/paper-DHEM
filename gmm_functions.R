library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)


## ------------------------------------------------------------
## Data generator + init
## ------------------------------------------------------------

invSigma_times <- function(R, v) {
  # v: length d
  # --- likely error point: dimension mismatch between v and R ---
  y <- backsolve(R, v, transpose = TRUE)  # solves R^T y = v
  backsolve(R, y, transpose = FALSE)      # solves R z = y  -> z = Sigma^{-1} v
}

rgmm <- function(n, pi, mu_list, Sigma_list) {
  K <- length(pi)
  d <- length(mu_list[[1]])
  z <- sample.int(K, size=n, replace=TRUE, prob=pi)
  x <- matrix(0, n, d)
  for (k in 1:K) {
    idx <- which(z == k)
    if (!length(idx)) next
    cholS <- chol(Sigma_list[[k]])
    Z <- matrix(rnorm(length(idx)*d), length(idx), d)
    x[idx, ] <- sweep(Z %*% cholS, 2, mu_list[[k]], "+")
  }
  list(x=x)
}

# Initialize GMM parameters: generate pi, mu, and Sigma with minimum separation constraint on means
gmm_init <- function(x, K,
                     pi_min = 0.05,
                     delta = 1,
                     Sigma_scale = 3,  
                     max_tries = 5000L,
                     seed = NULL
){
  if (!is.null(seed)) set.seed(seed)
  
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  
  if (K < 1) stop("K must be >= 1")
  if (pi_min < 0) stop("pi_min must be >= 0")
  if (K * pi_min > 1 + 1e-12) stop("Infeasible: K * pi_min must be <= 1")
  if (!is.finite(delta) || delta <= 0) stop("delta must be > 0")
  if (!is.finite(Sigma_scale) || Sigma_scale <= 0) stop("Sigma_scale must be > 0")
  
  # Sample covariance matrix
  S <- tryCatch(cov(x), error = function(e) diag(d))
  
  ## 1) pi: original scheme (ensure minimum pi_min for each component,
  ##    distribute the remaining mass randomly)
  if (K == 1) {
    pi <- 1
  } else {
    rem <- 1 - K * pi_min
    w <- rgamma(K, shape = 1, rate = 1)
    w <- w / sum(w)
    pi <- pi_min + rem * w
  }
  
  ## 2) Sigma: replicate the sample covariance matrix S for all K components
  Sigma_list <- replicate(K, S, simplify = FALSE)
  
  ## 3) mu: generate K points from N(mean(x), Sigma_scale * S),
  ##    repeat until all pairwise distances are at least delta
  xbar <- colMeans(x)
  Sig_mu <- Sigma_scale * S
  
  rmvnorm_one <- function(mu, Sigma) {
    # Base R only: sampling via Cholesky decomposition
    R <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(R)) {
      # If not numerically PSD, add diagonal jitter
      eps <- 1e-8
      Sigma2 <- Sigma + diag(eps, d)
      R <- chol(Sigma2)
    }
    as.numeric(mu + drop(t(R) %*% rnorm(d)))
  }
  
  mu_list <- vector("list", K)
  mu_list[[1]] <- as.numeric(xbar)
  
  for (k in 2:K) {
    ok <- FALSE
    for (t in 1:max_tries) {
      cand <- rmvnorm_one(xbar, Sig_mu)
      # Euclidean distance criterion: at least delta away from all previous mu's
      if (all(sapply(1:(k - 1), function(j) sum((cand - mu_list[[j]])^2) >= delta/10))) {
        mu_list[[k]] <- cand
        ok <- TRUE
        break
      }
    }
    if (!ok) stop(sprintf("Failed to generate mu with pairwise distance >= delta after %d tries (k=%d).",
                          max_tries, k))
  }
  
  list(pi = pi, mu = mu_list, Sigma = Sigma_list)
}


make_data_and_init <- function(n, K, theta_true, delta,seed=sample.int(1e9, 1)) {
  set.seed(seed)
  # cat("--- ",theta_true)
  x <- rgmm(n, theta_true$pi, theta_true$mu, theta_true$Sigma)$x
  
  set.seed(seed)
  theta_init <- gmm_init(x, K, delta = delta * 1.1)
  
  list(x = x, theta_init = theta_init,seed=seed)
}


dmvnorm_log <- function(X, m, S) {
  # X: n x d matrix
  # m: length d vector
  # S: d x d SPD matrix
  # return: length n vector of log N(x_i | m, S)
  
  X <- as.matrix(X)
  m <- as.numeric(m)
  d <- ncol(X)
  
  R <- chol(S)
  Xm <- sweep(X, 2, m, "-")                 # n x d
  Y  <- backsolve(R, t(Xm), transpose = TRUE) # d x n, solves R^T y = (x-m)^T
  quad <- colSums(Y^2)                      # length n
  logdet <- 2 * sum(log(diag(R)))
  -0.5 * (d * log(2 * pi) + logdet + quad)
}

gmm_Qr <- function(x, pi, mu, sigma, latent) {
  x <- as.matrix(x)
  n <- nrow(x)
  K <- length(pi)
  
  Q <- 0
  for (k in 1:K) {
    # log N(x_i | mu_k, Sigma_k) :n vector
    log_comp <- dmvnorm_log(x, mu[[k]], sigma[[k]])
    Q <- Q + sum(latent[, k] * (log(pi[k]) + log_comp))
  }
  Q
}



sumlog_mahalanobis_barrier_k <- function(k, mu, sigma, delta) {
  # b_k(mu_k) = log( sum_{l!=k} (mu_k-mu_l)' Sigma_k^{-1} (mu_k-mu_l) - delta )
  K <- length(mu)
  mu_k <- as.numeric(mu[[k]])
  Sig_k <- sigma[[k]]
  
  Rk <- chol(Sig_k)
  s <- 0
  
  for (l in 1:K) {
    if (l == k) next
    diff <- as.numeric(mu_k - mu[[l]])
    yk <- backsolve(Rk, diff, transpose = TRUE)  # Rk^T y = diff
    d2 <- sum(yk^2)                              # diff' Sig_k^{-1} diff
    s <- s + d2
  }
  
  inside <- s - delta
  if (inside <= 0) return(-Inf)
  log(inside)
}


gmm_BQr<- function(x, pi, mu, sigma, latent, bw, delta = 0.1) {
  x <- as.matrix(x)
  K <- length(pi)
  
  BQ <- 0
  for (k in 1:K) {
    log_comp <- dmvnorm_log(x, mu[[k]], sigma[[k]])
    BQ <- BQ + sum(latent[, k] * (log(pi[k]) + log_comp))
    
    bk <- sumlog_mahalanobis_barrier_k(k = k, mu = mu, sigma = sigma, delta = delta)
    if (!is.finite(bk)) return(-Inf)
    BQ <- BQ + bw * bk
  }
  BQ
}


solve_BQ_sumlog_mubarrier_newton <- function(k, x, mu, sigma, latent,
                                             delta, bw,
                                             maxiter = 100, tol = 1e-6,
                                             alpha_shrink = 0.5
) {
  x <- as.matrix(x)
  d <- ncol(x)
  K <- length(mu)
  
  # weights
  w <- latent[, k]
  Nk <- sum(w)
  if (Nk <= 0) return(as.numeric(mu[[k]]))
  
  s <- colSums(x * w)  # s_k
  
  # fixed Sigma_k and its chol
  Sig_k <- sigma[[k]]
  Rk <- chol(Sig_k)
  
  A_times <- function(v) invSigma_times(Rk, v)  # Sigma_k^{-1} v
  
  # other means fixed at this inner solve
  mu_others <- mu[-k]
  m <- Reduce(`+`, lapply(mu_others, as.numeric))
  Km1 <- K - 1L
  
  u_fun <- function(mk) Km1 * mk - m  # sum_{l!=k} (mk - mu_l)
  
  D_fun <- function(mk) {
    # D(mk) = sum_{l!=k} (mk-mu_l)' Sig_k^{-1} (mk-mu_l) - delta
    acc <- 0
    for (ml in mu_others) {
      diff <- mk - as.numeric(ml)
      Ad <- A_times(diff)
      acc <- acc + sum(diff * Ad)
    }
    acc - delta
  }
  
  G_fun <- function(mk) {
    u <- u_fun(mk)
    D <- D_fun(mk)
    D * (s - Nk * mk) + 2 * bw * u
  }
  
  J_fun <- function(mk) {
    # J_G(mk) = (-Nk*D + 2*bw*(K-1)) I + 2 (s - Nk*mk) (A u)^T
    u <- u_fun(mk)
    D <- D_fun(mk)
    Au <- A_times(u)
    
    base <- (-Nk * D + 2 * bw * Km1)
    base * diag(d) + 2 * ((s - Nk * mk) %*% t(Au))
  }
  
  mk <- as.numeric(mu[[k]])
  
  # Ensure initial feasibility: if not feasible, do small damping toward weighted mean
  if (D_fun(mk) <= 0) {
    mk_em <- s / Nk
    # try convex combination toward mk_em until feasible
    alpha <- 1
    for (bt in 1:30) {
      mk_try <- (1 - alpha) * mk + alpha * mk_em
      if (D_fun(mk_try) > 0) { mk <- mk_try; break }
      alpha <- alpha * 0.5
    }
    if (D_fun(mk) <= 0) return(as.numeric(mu[[k]]))  # fail-safe: no move
  }
  
  # Damped Newton iterations
  for (it in 1:maxiter) {
    G <- G_fun(mk)
    if (!all(is.finite(G))) break
    if (sqrt(sum(G^2)) < tol) break
    
    J <- J_fun(mk)
    if (!all(is.finite(J))) break
    
    step <- tryCatch(solve(J, -G), error = function(e) rep(0, d))
    
    # backtracking: keep D>0 and reduce ||G||
    alpha <- 1
    normG <- sqrt(sum(G^2))
    repeat {
      mk_try <- mk + alpha * step
      D_try <- D_fun(mk_try)
      if (D_try > 0) {
        G_try <- G_fun(mk_try)
        if (all(is.finite(G_try)) && sqrt(sum(G_try^2)) < normG) break
      }
      alpha <- alpha * alpha_shrink
      if (alpha < 1e-10) { mk_try <- mk; break }
    }
    
    mk <- mk_try
  }
  
  mk
}


is_feasible_sym <- function(mu, sigma, delta) {
  K <- length(mu)
  
  for (k in 1:(K - 1)) for (l in (k + 1):K) {
    diff <- as.numeric(mu[[k]] - mu[[l]])
    
    # --- likely error point: sigma[[k]] not SPD -> chol fails ---
    Rk <- chol(sigma[[k]])
    yk <- backsolve(Rk, diff, transpose = TRUE)
    d2k <- sum(yk^2)
    
    # --- likely error point: sigma[[l]] not SPD -> chol fails ---
    Rl <- chol(sigma[[l]])
    yl <- backsolve(Rl, diff, transpose = TRUE)
    d2l <- sum(yl^2)
    
    if ((d2k + d2l) <= delta) return(FALSE)
  }
  TRUE
}

is_feasible_sumlog <- function(mu, sigma, delta) {
  K <- length(mu)
  for (k in 1:K) {
    val <- sumlog_mahalanobis_barrier_k(k, mu, sigma, delta)
    if (!is.finite(val)) return(FALSE)
  }
  TRUE
}

# Accept a proposed (mu, Sigma) update only if it is feasible and does not decrease the barrier-augmented objective
accept_mu_sigma_by_damping_feasible <- function(mu_old, sigma_old,
                                                x, pi_new, latent,
                                                mu_prop,
                                                bw, delta,
                                                BQ_old,
                                                alpha0 = 1.0,
                                                alpha_shrink = 0.5,
                                                max_backtrack = 20
) {
  K <- length(mu_old)
  
  # compute Sigma and objective for the proposal
  sigma_prop <- gmm_sigma_update(x, mu_prop, latent)
  
  # check feasibility of the proposal
  feasible_prop <- is_feasible_sym(mu_prop, sigma_prop, delta)
  
  if (feasible_prop) {
    BQ_prop <- gmm_BQr(x, pi_new, mu_prop, sigma_prop,
                       latent = latent, bw = bw, delta = delta)
    
    # accept immediately if feasible and non-decreasing
    if (is.finite(BQ_old) && is.finite(BQ_prop) && (BQ_prop >= BQ_old)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
    
    # accept if the old objective is not finite but the proposal is finite
    if (!is.finite(BQ_old) && is.finite(BQ_prop)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
  }
  
  # backtracking: interpolate between old and proposed mu
  alpha <- alpha0
  mu_try <- mu_old
  sigma_try <- sigma_old
  BQ_try <- BQ_old
  
  for (bt in 1:max_backtrack) {
    alpha <- alpha * alpha_shrink
    
    mu_try <- lapply(seq_len(K), function(k) {
      (1 - alpha) * as.numeric(mu_old[[k]]) + alpha * as.numeric(mu_prop[[k]])
    })
    
    sigma_try <- gmm_sigma_update(x, mu_try, latent)
    
    # skip infeasible candidates
    if (!is_feasible_sym(mu_try, sigma_try, delta)) next
    
    BQ_try <- gmm_BQr(x, pi_new, mu_try, sigma_try,
                      latent = latent, bw = bw, delta = delta)
    
    if (is.finite(BQ_try) && is.finite(BQ_old) && (BQ_try >= BQ_old)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
    
    if (!is.finite(BQ_old) && is.finite(BQ_try)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
  }
  
  # fallback: keep the old parameters
  return(list(mu = mu_old, sigma = sigma_old, BQ = BQ_old, alpha = 0.0))
}


gmm_B_value <- function(theta, delta = 0.1, symmetric = TRUE) {
  mu    <- theta$mu
  sigma <- theta$Sigma
  K <- length(mu)
  
  B <- 0
  for (k in seq_len(K)) {
    bk = sumlog_mahalanobis_barrier_k(k = k, mu = mu, sigma = sigma, delta = delta)
    B <- B + bk
  }
  B
}

# Delta |B(theta0, theta1)| = | B(theta1) - B(theta0) |
gmm_delta_abs_B <- function(theta0, theta1, delta = 0.1, symmetric = TRUE) {
  B0 <- gmm_B_value(theta0, delta = delta, symmetric = symmetric)
  B1 <- gmm_B_value(theta1, delta = delta, symmetric = symmetric)
  
  if (!is.finite(B0) || !is.finite(B1)) return(Inf)
  
  abs(B1 - B0)
}


# Cap the barrier weight so that the barrier change does not exceed the KL-based target
cap_bw_by_barrier <- function(theta0, theta1, eta, dkl_1, bw, delta, symmetric = TRUE) {
  absDeltaB <- gmm_delta_abs_B(theta0, theta1, delta = delta, symmetric = symmetric)
  
  # no restriction if the barrier change is zero or not finite
  if (!is.finite(absDeltaB) || absDeltaB <= 0) {
    return(list(bw_new = bw, violated = FALSE, absDeltaB = absDeltaB))
  }
  
  delta_target <- eta * dkl_1
  
  bw_cap <- delta_target / absDeltaB
  bw_new <- min(bw, bw_cap)
  
  violated <- (bw > bw_cap)
  list(bw_new = bw_new, violated = violated, absDeltaB = absDeltaB, bw_cap = bw_cap, delta_target = delta_target)
}
# ------------------------------------------------------------------
# Sumlog-barrier version of accept_mu_sigma_by_damping_feasible()
#   - feasibility: is_feasible_sumlog()
#   - objective  : gmm_BQr()
# ------------------------------------------------------------------
accept_mu_sigma_by_damping_feasible_sumlog <- function(mu_old, sigma_old,
                                                       x, pi_new, latent,
                                                       mu_prop,
                                                       bw, delta,
                                                       BQ_old,
                                                       alpha0 = 1.0,
                                                       alpha_shrink = 0.5,
                                                       max_backtrack = 20) {
  K <- length(mu_old)
  
  # 1) proposal sigma + feasibility + BQ
  sigma_prop <- gmm_sigma_update(x, mu_prop, latent)
  
  # sumlog feasibility (all k: sum_{l!=k} d2_{kl}(Sigma_k^{-1}) - delta > 0)
  feasible_prop <- is_feasible_sumlog(mu_prop, sigma_prop, delta)
  
  if (feasible_prop) {
    BQ_prop <- gmm_BQr(x, pi_new, mu_prop, sigma_prop,
                       latent = latent, bw = bw, delta = delta)
    
    # feasible + non-decrease -> accept
    if (is.finite(BQ_old) && is.finite(BQ_prop) && (BQ_prop >= BQ_old)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
    # if old is not finite but new is finite, accept
    if (!is.finite(BQ_old) && is.finite(BQ_prop)) {
      return(list(mu = mu_prop, sigma = sigma_prop, BQ = BQ_prop, alpha = 1.0))
    }
  }
  
  # 2) backtracking: mu_try = (1-a)mu_old + a*mu_prop
  alpha <- alpha0
  mu_try <- mu_old
  sigma_try <- sigma_old
  BQ_try <- BQ_old
  
  for (bt in 1:max_backtrack) {
    alpha <- alpha * alpha_shrink
    
    mu_try <- lapply(seq_len(K), function(k) {
      (1 - alpha) * as.numeric(mu_old[[k]]) + alpha * as.numeric(mu_prop[[k]])
    })
    
    sigma_try <- gmm_sigma_update(x, mu_try, latent)
    
    # feasibility first
    if (!is_feasible_sumlog(mu_try, sigma_try, delta)) next
    
    BQ_try <- gmm_BQr(x, pi_new, mu_try, sigma_try,
                      latent = latent, bw = bw, delta = delta)
    
    if (is.finite(BQ_try) && is.finite(BQ_old) && (BQ_try >= BQ_old)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
    if (!is.finite(BQ_old) && is.finite(BQ_try)) {
      return(list(mu = mu_try, sigma = sigma_try, BQ = BQ_try, alpha = alpha))
    }
  }
  
  # 3) fallback: null move
  return(list(mu = mu_old, sigma = sigma_old, BQ = BQ_old, alpha = 0.0))
}


## ------------------------------------------------------------
## Helper func
## ------------------------------------------------------------

gmm_estep_annealed <- function(x, theta, r = 1) {
  x <- as.matrix(x)
  n <- nrow(x)
  K <- length(theta$pi)
  
  # log_resp[i,k] = r * (log pi_k + log N(x_i | mu_k, Sigma_k))
  log_resp <- matrix(NA_real_, nrow = n, ncol = K)
  
  for (k in 1:K) {
    log_comp <- dmvnorm_log(x, theta$mu[[k]], theta$Sigma[[k]])
    log_resp[, k] <- r * (log(theta$pi[k]) + log_comp)
  }
  
  # w_ik = exp(log_resp_ik)
  w <- exp(log_resp)
  
  # gamma_ik = w_ik / sum_j w_ij
  gamma <- w / rowSums(w)
  
  gamma
}

gmm_pi_update <- function(x, latent){
  N <- nrow(x)
  pi <- colSums(latent) / N
  as.numeric(pi)
}

gmm_mu_update <- function(x, latent){
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  K <- ncol(latent)
  
  Nk <- colSums(latent)          # K-vector
  mu <- vector("list", K)        
  
  for(k in 1:K){
    w <- latent[, k]             # n-vector
    mu[[k]] <- as.numeric(colSums(x * w) / Nk[k])
  }
  
  mu
}

gmm_sigma_update <- function(x, mu, latent){
  x <- as.matrix(x)
  d <- ncol(x)
  K <- ncol(latent)
  
  Nk <- colSums(latent)
  Sigma <- vector("list", K)
  
  for(k in 1:K){
    xc <- sweep(x, 2, mu[[k]], "-")     # n x d
    w  <- latent[, k]                   # n-vector
    Sigma[[k]] <- (t(xc) %*% (xc * w)) / Nk[k]  # d x d
  }
  
  Sigma
}

gmm_em_at_r <- function(x, theta_init, r = 1,
                        tol = 1e-3, max_iter = 200, verbose = FALSE) {
  theta <- theta_init
  trace <- data.frame(iter=integer(0), r=numeric(0), Qr=numeric(0))
  
  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma
  
  Qr_old <- gmm_Qr(x, pi, mu, sigma, latent=gmm_estep_annealed(x, theta, r))
  
  # write trace
  trace <- rbind(trace, data.frame(iter=0, r=r, Qr=Qr_old))
  
  for (it in 1:max_iter) {
    # E-step
    latent = gmm_estep_annealed(x, theta, r)
    
    # M-step
    pi_new    = gmm_pi_update(x,latent)
    mu_new    = gmm_mu_update(x,latent)
    sigma_new = gmm_sigma_update(x,mu_new,latent)
    
    # Comupte new Qr  
    Qr_new <- gmm_Qr(x, pi_new, mu_new, sigma_new, latent)
    
    # Update parameter
    theta$pi <- pi_new
    theta$mu <- mu_new
    theta$Sigma <- sigma_new
    
    # write trace 
    trace <- rbind(trace, data.frame(iter=it, r=r, Qr=Qr_new))
    if (verbose) cat("iter:", it, "Qr:", Qr_new, "dQr:", Qr_new - Qr_old, "\n")
    
    # Break rule 
    if(abs(Qr_old-Qr_new)<tol){
      break}else{Qr_old=Qr_new}
  }
  
  
  list(theta = theta, trace = trace)
}

gmm_em_at_r_bw <- function(x, theta,
                           r, bw,
                           delta,
                           tol = 1e-8, max_iter = 1000,
                           verbose = FALSE) {
  # optimization settings for inner mu solver
  maxiter_optim <- 200
  tol_optim <- 1e-4
  
  # trace of objective values
  trace <- data.frame(iter=integer(0), r=numeric(0), bw=numeric(0), Qr=numeric(0))
  
  # unpack parameters
  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma
  
  K <- length(mu)
  mu_new <- vector("list", K)
  
  ## initial barrier-augmented Q
  BQr_old <- gmm_BQr(
    x, pi, mu, sigma,
    latent = gmm_estep_annealed(x, theta, r),
    bw = bw, delta = delta
  )
  
  trace <- rbind(trace, data.frame(iter=0, r=r, bw=bw, Qr=BQr_old))
  
  for (it in 1:max_iter) {
    ## E-step
    latent <- gmm_estep_annealed(x, theta, r)
    
    ## M-step: update pi
    pi_new <- gmm_pi_update(x, latent)
    
    ## M-step: update mu by Gauss–Seidel with sumlog barrier solver
    mu_prop <- mu
    for (k in seq_len(K)) {
      mu_prop[[k]] <- solve_BQ_sumlog_mubarrier_newton(
        k = k, x = x,
        mu = mu_prop,
        sigma = sigma,
        latent = latent,
        delta = delta, bw = bw,
        maxiter = maxiter_optim, tol = tol_optim,
        alpha_shrink = 0.5
      )
    }
    
    ## accept update with feasibility and BQ non-decrease check
    acc <- accept_mu_sigma_by_damping_feasible_sumlog(
      mu_old = mu, sigma_old = sigma,
      x = x, pi_new = pi_new, latent = latent,
      mu_prop = mu_prop,
      bw = bw, delta = delta,
      BQ_old = BQr_old,
      alpha0 = 1.0, alpha_shrink = 0.5, max_backtrack = 20
    )
    
    mu_new    <- acc$mu
    sigma_new <- acc$sigma
    BQr_new   <- acc$BQ
    
    ## update parameters
    theta$pi    <- pi_new
    theta$mu    <- mu_new
    theta$Sigma <- sigma_new
    
    ## refresh locals
    pi    <- pi_new
    mu    <- mu_new
    sigma <- sigma_new
    
    ## update trace
    trace <- rbind(trace, data.frame(iter=it, r=r, bw=bw, Qr=BQr_new))
    
    if (verbose) {
      cat("iter:", it,
          "BQr(sumlog):", BQr_new,
          "dBQr:", BQr_new - BQr_old,
          "\n")
    }
    
    ## stop if improvement is small
    if (abs(BQr_old - BQr_new) < tol) {
      break
    } else {
      BQr_old <- BQr_new
    }
  }
  
  list(theta = theta, trace = trace)
}

make_trace_entry <- function(step, out, trace_full = TRUE, ...) {
  meta <- list(...)
  if (trace_full) {
    return(c(list(step = step), meta, list(out = out)))
  } else {
    return(c(list(
      step = step,
      BQ_final = tail(out$trace$Qr, 1),
      iters = nrow(out$trace) - 1
    ), meta))
  }
}

gmm_delta_DKL <- function(x, theta0, theta1, r) {
  x <- as.matrix(x)
  n <- nrow(x)
  
  # standard posterior (r=1) for theta0, theta1
  gamma0 <- gmm_estep_annealed(x, theta0, r = 1)
  gamma1 <- gmm_estep_annealed(x, theta1, r = 1)
  
  # annealed posterior weight (r=r) under theta0
  w <- gmm_estep_annealed(x, theta0, r = r)
  
  # per-i contribution: sum_k w_ik * (log gamma0_ik - log gamma1_ik)
  contrib <- rowSums(w * (log(gamma0) - log(gamma1)))
  
  # empirical mean over i
  mean(contrib)
}


gmm_find_adap_dhem_at_r <- function(x, theta, r = 1,
                                    bw = 1e-2,
                                    eta = 0.1,
                                    delta=0.1,
                                    tol = 1e-3, max_iter = 200,
                                    verbose = FALSE) {
  
  # optimization settings for inner mu solver
  maxiter_optim <- 200
  tol_optim <- 1e-4
  
  # trace of objective values
  trace <- data.frame(iter=integer(0), r=numeric(0), bw=numeric(0), Qr=numeric(0))
  
  # store initial theta for fallback
  theta0 = theta
  
  # unpack parameters
  pi    <- theta$pi
  mu    <- theta$mu
  sigma <- theta$Sigma
  
  K <- length(mu)
  mu_new <- vector("list", K)
  
  ## initial BQ at theta0
  BQr_old <- gmm_BQr(
    x, pi, mu, sigma,
    latent = gmm_estep_annealed(x, theta, r),
    bw = bw, delta = delta
  )
  
  trace <- rbind(trace, data.frame(iter=0, r=r, bw=bw, Qr=BQr_old))
  
  for (it in 1:max_iter) {
    
    ## ----- E-step -----
    latent <- gmm_estep_annealed(x, theta, r)
    
    ## ----- M-step -----
    pi_new <- gmm_pi_update(x, latent)
    
    ## ----- M-step: mu update (Gauss–Seidel with sumlog barrier solver) -----
    mu_prop <- mu
    for (k in seq_len(K)) {
      mu_prop[[k]] <- solve_BQ_sumlog_mubarrier_newton(
        k = k, x = x,
        mu = mu_prop,
        sigma = sigma,
        latent = latent,
        delta = delta, bw = bw,
        maxiter = maxiter_optim, tol = tol_optim,
        alpha_shrink = 0.5
      )
    }
    
    ## ----- Backtracking acceptance (feasible + BQ non-decrease) -----
    acc <- accept_mu_sigma_by_damping_feasible_sumlog(
      mu_old = mu, sigma_old = sigma,
      x = x, pi_new = pi_new, latent = latent,
      mu_prop = mu_prop,
      bw = bw, delta = delta,
      BQ_old = BQr_old,
      alpha0 = 1.0, alpha_shrink = 0.5, max_backtrack = 20
    )
    
    mu_new    <- acc$mu
    sigma_new <- acc$sigma
    BQr_new   <- acc$BQ
    
    ## update parameters
    theta$pi    <- pi_new
    theta$mu    <- mu_new
    theta$Sigma <- sigma_new
    
    ## refresh locals
    pi    <- pi_new
    mu    <- mu_new
    sigma <- sigma_new
    
    ## break if infeasible
    feasible_now <- is.finite(BQr_new)
    if(!is.finite(BQr_new)){break}
    
    ## update trace
    trace <- rbind(trace, data.frame(iter=it, r=r, bw=bw, Qr=BQr_new))
    
    ## ----- Acceptance condition (adaptive criterion) -----
    dkl_r <- gmm_delta_DKL(x, theta0, theta, r = r)
    dkl_1 <- gmm_delta_DKL(x, theta0, theta, r = 1)
    target <- eta * dkl_1
    acc_now <- (dkl_r >= target)
    
    trace <- rbind(trace,data.frame(iter = it, r = r, bw = bw,Qr = BQr_new))
    
    if (verbose) {
      cat("iter:", it,
          "BQ:", BQr_new,
          "dBQ:", BQr_new - BQr_old,
          "feasible:", feasible_now,
          "dkl_r:", dkl_r,
          "target:", target,
          "accept:", acc_now, "\n")
    }
    
    ## return if acceptance condition satisfied
    if (acc_now) {
      return(list(theta = theta,
                  trace = trace,
                  accept = TRUE,
                  bw = bw,
                  BQ_old = BQr_old,
                  BQ_new = BQr_new,
                  dkl_r = dkl_r,
                  dkl_1 = dkl_1,
                  target = target))
    }
    
    ## stop if improvement is small
    if (abs(BQr_new - BQr_old) < tol) break
    
    ## update for next iteration
    BQr_old <- BQr_new
  }
  
  ## return initial theta if no acceptable update found
  list(theta = theta0,
       trace = trace,
       accept = FALSE,
       bw = bw,
       BQ_old = BQr_old,
       BQ_new = BQr_new,
       dkl_r = dkl_r,
       dkl_1 = dkl_1,
       target = target)
}




## ------------------------------------------------------------
## Algorithms
## ------------------------------------------------------------
gmm_em <- function(x, theta_init,
                   tol = 1e-4, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  out = gmm_em_at_r(x, theta_init, r = 1, tol = tol, max_iter = max_iter, verbose = verbose)
  
  trace <- make_trace_entry(step = 1, out = out, trace_full = trace_full, r = 1)
  
  list(theta = out$theta, trace = trace)
}

gmm_daem <- function(x, theta_init,
                     r_init = 0.2, n_steps = 50,
                     tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  theta <- theta_init
  trace <- vector("list", n_steps)
  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))
  
  for (s in seq_along(r_grid)) {
    r <- r_grid[s]
    out <- gmm_em_at_r(x, theta, r = r, tol = tol, max_iter = max_iter, verbose = verbose)
    theta <- out$theta
    
    # r by write
    trace[[s]] <- make_trace_entry(step = s, out = out, trace_full = trace_full, r = r)
  }
  list( theta = theta,trace = trace)
}


gmm_dhem <- function(x, theta_init,
                     r_init = 0.2, n_steps = 50,
                     bw_init = 1, bw_end = 1e-10,
                     delta = 0.1,
                     tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  
  theta <- theta_init
  trace <- vector("list", n_steps)
  
  r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))                 # r -> 1
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))      # bw -> bw_end
  
  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]
    bw <- bw_grid[s]
    
    out <- gmm_em_at_r_bw(
      x, theta,
      r = r, bw = bw, delta = delta,
      tol = tol, max_iter = max_iter, verbose = verbose
    )
    
    theta <- out$theta
    
    trace[[s]] <- make_trace_entry(step = s, out = out, trace_full = trace_full, r = r, bw = bw)
  }
  
  list(theta = theta, trace = trace  )
}

gmm_barrier <- function(x, theta_init,
                        n_steps = 50,
                        # r_init = 0.2,
                        bw_init = 1, bw_end = 1e-10,
                        delta = 0.1,
                        tol = 1e-8, max_iter = 200, verbose = FALSE,trace_full = TRUE) {
  
  theta <- theta_init
  trace <- vector("list", n_steps)
  
  # r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))               
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))      
  
  r = 1
  for (s in seq_len(n_steps)) {
    # r  <- r_grid[s]
    bw <- bw_grid[s]
    
    out <- gmm_em_at_r_bw(
      x, theta,
      r = r, bw = bw, delta = delta,
      tol = tol, max_iter = max_iter, verbose = verbose
    )
    
    theta <- out$theta
    
    trace[[s]] <- make_trace_entry(step=s, out=out, trace_full=trace_full, r=1, bw=bw)
  }
  
  list(theta = theta, trace = trace  )
}



gmm_adaptive_dhem <- function(x, theta,
                              r_init = 0.2, n_steps = 50,
                              bw_init = 1e-2,
                              delta = 0.1,
                              eta = 0.1,
                              tol_BQ = 1e-4,
                              tol_find = 1e-3,
                              max_iter_inner = 50,
                              max_iter_find = 200,
                              verbose = FALSE,
                              trace_full = TRUE) {
  
  x <- as.matrix(x)
  
  # geometric schedule for r
  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))
  
  # barrier weight (monotonically decreasing)
  bw <- bw_init
  
  trace <- vector("list", n_steps)
  
  for (s in seq_along(r_grid)) {
    
    r <- r_grid[s]
    
    if (verbose) cat("\n===== step", s, " r =", r, " bw =", bw, "=====\n")
    
    BQ_prev_acc <- -Inf
    out_last <- NULL
    acc_count <- 0L
    
    # inner loop at fixed r
    for (inner in seq_len(max_iter_inner)) {
      
      out <- gmm_find_adap_dhem_at_r(
        x = x,
        theta = theta,
        r = r,
        bw = bw,
        eta = eta,
        delta = delta,
        tol = tol_find,
        max_iter = max_iter_find,
        verbose = verbose
      )
      out_last <- out
      
      # accept1 failure → move to next r
      if (!out$accept) {
        if (verbose) cat("  inner:", inner, " accept1=FALSE -> next r\n")
        break
      }
      
      theta_new <- out$theta
      BQ_new    <- out$BQ_new
      
      # check barrier condition (2)
      cap <- cap_bw_by_barrier(
        theta0 = theta,
        theta1 = theta_new,
        eta = eta,
        dkl_1 = out$dkl_1,
        bw = bw,
        delta = delta
      )
      
      # shrink bw if violated and retry at same r
      if (isTRUE(cap$violated)) {
        if (verbose) {
          cat("  inner:", inner,
              " barrier violated -> bw:",
              bw, "->", cap$bw_new, "\n")
        }
        bw <- cap$bw_new
        next
      }
      
      # accept update
      theta <- theta_new
      acc_count <- acc_count + 1L
      
      # stop if BQ improvement is small
      if (is.finite(BQ_new) && abs(BQ_new - BQ_prev_acc) < tol_BQ) {
        if (verbose) cat("  inner:", inner, " small ΔBQ -> break\n")
        BQ_prev_acc <- BQ_new
        break
      }
      BQ_prev_acc <- BQ_new
    }
    
    # store trace
    if (trace_full) {
      trace[[s]] <- list(step = s, r = r, bw = bw,
                         accepted_moves = acc_count,
                         out_last = out_last,
                         theta = theta)
    } else {
      if (is.null(out_last)) {
        trace[[s]] <- list(step = s, r = r, bw = bw,
                           accepted_moves = 0L,
                           accept_last = FALSE,
                           BQ_final = NA_real_,
                           iters_last = 0L)
      } else {
        BQ_final <- tail(out_last$trace$BQ, 1)
        trace[[s]] <- list(step = s, r = r, bw = bw,
                           accepted_moves = acc_count,
                           accept_last = out_last$accept,
                           BQ_final = BQ_final,
                           iters_last = nrow(out_last$trace) - 1L)
      }
    }
  }
  
  list(theta = theta, bw_final = bw, trace = trace)
}



## -----------------------------
## distances + matching (minimal)
## -----------------------------
l2_dist <- function(a, b) sqrt(sum((as.numeric(a) - as.numeric(b))^2))
fro_dist <- function(A, B) sqrt(sum((A - B)^2))

match_by_mu <- function(theta_hat, theta_true) {
  K <- length(theta_true$pi)
  C <- matrix(0, K, K)
  for (k in 1:K) for (j in 1:K) {
    C[k, j] <- l2_dist(theta_hat$mu[[k]], theta_true$mu[[j]])
  }
  used_hat <- rep(FALSE, K)
  perm_hat_to_true <- integer(K)
  for (j in 1:K) {
    cand <- which(!used_hat)
    kbest <- cand[which.min(C[cand, j])]
    perm_hat_to_true[j] <- kbest
    used_hat[kbest] <- TRUE
  }
  th <- list(
    pi = theta_hat$pi[perm_hat_to_true],
    mu = theta_hat$mu[perm_hat_to_true],
    Sigma = theta_hat$Sigma[perm_hat_to_true]
  )
  list(theta = th, perm = perm_hat_to_true)
}

compute_err3 <- function(theta_hat, theta_true) {
  K <- length(theta_true$pi)
  th <- match_by_mu(theta_hat, theta_true)$theta
  
  pi_abs <- mean(abs(th$pi - theta_true$pi))
  
  mu_L2 <- mean(vapply(seq_len(K),
                       function(j) l2_dist(th$mu[[j]], theta_true$mu[[j]]),
                       numeric(1)))
  
  Sigma_F <- mean(vapply(seq_len(K),
                         function(j) fro_dist(th$Sigma[[j]], theta_true$Sigma[[j]]),
                         numeric(1)))
  
  c(pi = pi_abs, mu = mu_L2, Sigma = Sigma_F)
}

## -----------------------------
## one simulation summary (requested format)
## -----------------------------

one_sim_summary <- function(simulation_num = 1,
                            n, K, theta_true,
                            tol_limit = 1e-6,
                            r_init = 0.2, n_steps = 20,
                            bw_init = 1e-2, bw_end = 1e-10,
                            delta = 0.5,
                            max_iter = 1000,eta=0.1,
                            verbose = FALSE) {
  
  # random seed for this simulation
  seed <- sample.int(1e9, 1)
  
  # generate data and initial parameters
  tmp <- make_data_and_init(n, K, theta_true, delta)
  x <- tmp$x
  theta0 <- tmp$theta_init
  
  # safe wrapper: return NULL if error occurs
  safe_run <- function(expr) tryCatch(expr, error = function(e) NULL)
  
  # run all methods
  res <- list(
    EM = safe_run(gmm_em(x, theta0, tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    DAEM = safe_run(gmm_daem(x, theta0,
                             r_init = r_init, n_steps = n_steps,
                             tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    DHEM = safe_run(gmm_dhem(x, theta0,
                             r_init = r_init, n_steps = n_steps,
                             bw_init = bw_init, bw_end = bw_end,
                             delta = delta,
                             tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    Barrier = safe_run(gmm_barrier(x, theta0,
                                   n_steps = n_steps,
                                   bw_init = bw_init, bw_end = bw_end,
                                   delta = delta,
                                   tol = tol_limit, max_iter = max_iter, verbose = verbose)),
    AdaptiveDHEM = safe_run(gmm_adaptive_dhem(x, theta0,
                                              r_init = r_init, n_steps = n_steps,
                                              bw_init = bw_init,
                                              delta = delta,
                                              eta = 0.1,
                                              tol_BQ = tol_limit, tol_find = 1e-3,
                                              max_iter_inner = max_iter, max_iter_find = 200,
                                              verbose = verbose, trace_full = TRUE))
  )
  
  # extract estimated parameters
  theta_list <- lapply(res, function(o) if (is.null(o)) NULL else o$theta)
  
  methods <- names(theta_list)
  paras <- c("pi", "mu", "Sigma")
  
  rows <- list()
  idx <- 0L
  
  for (m in methods) {
    th_hat <- theta_list[[m]]
    succ <- if (is.null(th_hat)) 0L else 1L
    
    # if failed, record NA errors
    if (succ == 0L) {
      for (p in paras) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          simulation_num = simulation_num,
          seed = seed,
          method = m,
          para = p,
          err = NA_real_,
          succ = succ,
          stringsAsFactors = FALSE
        )
      }
    } else {
      # compute estimation error for each parameter
      e <- compute_err3(th_hat, theta_true)
      for (p in paras) {
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          simulation_num = simulation_num,
          seed = seed,
          method = m,
          para = p,
          err = as.numeric(e[p]),
          succ = succ,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # combine results into a single data frame
  do.call(rbind, rows)
}
