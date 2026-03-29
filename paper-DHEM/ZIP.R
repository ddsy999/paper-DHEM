# ZIP data generator ------------------------------------------------------

rzip <- function(n, pi, lambda) {
  # Z=0 (structural zero) with prob pi
  # Z=1 (Poisson) with prob 1-pi
  z <- rbinom(n, size = 1, prob = 1 - pi)  # z=1 => Poisson component
  x <- z * rpois(n, lambda = lambda)       # if z=0 => 0
  x
}

compute_lambda_init <- function(y, set_pi_init) {
  w <- ifelse(y == 0, set_pi_init, 1)
  sum(w * y) / sum(w)
}

rel_err <- function(est, true) {
  abs(est - true) / true
}

abs_err <- function(est, true) {
  abs(est - true) 
}

err <- function(est, true) {
  (est - true) 
}

ratio_err <- function(est, true) {
  max(abs(est/true),1/abs(est/true))
}

log_err <- function(est, true) {
  abs(log(est)-log(true))
}


# Annealed ZIP E-step ------------------------------------------------------

zip_estep_annealed <- function(x, pi, lambda, r=1) {
  if (r <= 0) stop("r must be > 0")
  
  gamma <- numeric(length(x))
  
  # x > 0 : 반드시 Poisson 성분
  pos <- (x > 0L)
  gamma[pos] <- 1.0
  
  # x = 0 : annealed posterior
  z0 <- !pos
  if (any(z0)) {
    w1 <- (1 - pi) * exp(-lambda)  # Poisson에서 0이 나올 joint weight
    w0 <- pi                       # structural zero weight
    
    w1r <- w1^r
    w0r <- w0^r
    gamma[z0] <- w1r / (w0r + w1r)
  }
  # γ=(1−π)e−λ/(π+(1−π)e−λ)
  gamma
}

zip_mstep <- function(x, gamma) {
  n <- length(x)
  A <- sum(1 - gamma)  # structural-zero 기대 개수
  B <- sum(gamma)      # Poisson component 기대 개수
  
  pi_new <- A / n
  lambda_new <- sum(gamma * x) / B
  
  list(pi = pi_new, lambda = lambda_new)
}


zip_em <- function(x,
                   pi_init = 0.5,
                   lambda_init = NULL,
                   tol = 1e-10,
                   max_iter = 200
                  ) {
  n <- length(x)
  
  # 초기값
  pi <- pi_init
  if (is.null(lambda_init)) {
    lambda <- compute_lambda_init(x, set_pi_init = pi_init)
  } else {
    lambda <- lambda_init
  }
  
  trace <- data.frame(
    iter = integer(0),
    pi = numeric(0),
    lambda = numeric(0),
    d_pi = numeric(0),
    d_lambda = numeric(0)
  )
  
  for (it in seq_len(max_iter)) {
    # E-step (r=1)
    gamma <- zip_estep_annealed(x, pi = pi, lambda = lambda, r = 1)
    
    # M-step (사용자 제공 함수 그대로)
    th_new <- zip_mstep(x, gamma)
    pi_new <- th_new$pi
    lambda_new <- th_new$lambda
    
    # 변화량 기록
    d_pi <- abs(pi_new - pi)
    d_lam <- abs(lambda_new - lambda)
    
    trace[it, ] <- list(it, pi_new, lambda_new, d_pi, d_lam)
    
    # 수렴 체크
    if (max(d_pi, d_lam) < tol) {
      pi <- pi_new
      lambda <- lambda_new
      break
    }
    
    # 업데이트
    pi <- pi_new
    lambda <- lambda_new
  }
  
  list(
    theta = list(pi = pi, lambda = lambda),
    trace = trace
  )
}
# ============================================================
# 0) 기존 함수들이 이미 정의되어 있다고 가정:
#    - zip_estep_annealed(x, pi, lambda, r)
#    - zip_mstep(x, gamma)
#    - zip_dkl_std_post(x, theta0, theta1)
#    - zip_delta_dkl(x, theta0, theta1, r)
# ============================================================


# ============================================================
# 1) fixed-r에서 EM을 tol까지 돌리는 "wrapper" (새 유틸)
#    -> 기존 zip_estep_annealed + zip_mstep만 사용
# ============================================================
zip_em_at_r <- function(x, theta_init, r,
                        tol = 1e-10,
                        max_iter = 200) {
  theta <- theta_init
  
  trace <- data.frame(
    iter   = integer(0),
    r      = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (it in 1:max_iter) {
    gamma  <- zip_estep_annealed(x, pi = theta$pi, lambda = theta$lambda, r = r)
    theta1 <- zip_mstep(x, gamma)
    
    d_pi  <- abs(theta1$pi - theta$pi)
    d_lam <- abs(theta1$lambda - theta$lambda)
    
    theta <- theta1
    
    trace <- rbind(trace, data.frame(
      iter   = it,
      r      = r,
      pi     = theta$pi,
      lambda = theta$lambda,
      d_pi   = d_pi,
      d_lam  = d_lam
    ))
    
    if (max(d_pi, d_lam) < tol) break
  }
  
  list(theta = theta, trace = trace)
}


# ============================================================
# 2) DAEM: r-grid (log schedule) 각 r에서 zip_em_at_r로 수렴
#    -> 기존 zip_daem의 "1회 업데이트" 문제 해결
# ============================================================
zip_daem <- function(x,
                     theta_init,
                     r_init  = 0.2,
                     n_steps = 50,
                     tol = 1e-10,
                     max_iter = 200) {
  
  r_grid <- exp(seq(log(r_init), 0, length.out = n_steps))  # ends at 1
  
  theta <- theta_init
  
  trace <- data.frame(
    step   = integer(0),
    r      = numeric(0),
    iter   = integer(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (s in 1:n_steps) {
    r <- r_grid[s]
    
    res  <- zip_em_at_r(x, theta_init = theta, r = r, tol = tol, max_iter = max_iter)
    theta <- res$theta
    
    tr <- res$trace
    tr$step <- s
    trace <- rbind(trace, tr[, c("step","r","iter","pi","lambda","d_pi","d_lam")])
  }
  
  list(theta = theta, trace = trace, r_grid = r_grid)
}


# ============================================================
# 3) Adaptive DAEM:
#    - r <= r_switch: fixed phase (각 r에서 inner EM 수렴, 항상 accept)
#    - r >  r_switch: 후보 r들을 증가시키며 (rescan)
#         theta1(r) = inner EM 수렴 결과로 만들고
#         zip_delta_dkl >= eta * zip_dkl_std_post 이면 accept
#         아니면 다음 r로 넘어감
# ============================================================

zip_BQ <- function(x, theta, gamma, bw = 0, pi_min = 0) {
  pi <- theta$pi
  lam <- theta$lambda
  if (!is.finite(pi) || !is.finite(lam)) return(-Inf)
  if (pi <= 0 || pi >= 1 || lam <= 0) return(-Inf)
  if (bw > 0 && pi <= pi_min) return(-Inf)
  
  Q <- sum((1 - gamma) * log(pi) + gamma * log(1 - pi) +
             gamma * (x * log(lam) - lam))
  B <- if (bw > 0) bw * log(pi - pi_min) else 0
  Q + B
}


# --- Barrier pi-update (M-step subroutine) -------------------------------

zip_mstep_pi_barrier <- function(A, B, bw, pi_min) {
  # maximize: f(pi) = A log(pi) + B log(1-pi) + bw log(pi - pi_min)
  # domain: pi in (pi_min, 1)
  # FOC -> quadratic:
  # n*pi^2 - (n*(1+pi_min)+bw)*pi + A*pi_min = 0
  n <- A + B
  
  b <- n * (1 + pi_min) + bw
  c <- A * pi_min
  
  disc <- b*b - 4*n*c
  if (disc < 0) disc <- 0  # numeric guard
  
  r1 <- (b - sqrt(disc)) / (2*n)
  r2 <- (b + sqrt(disc)) / (2*n)
  
  eps <- 1e-12
  lo  <- pi_min + eps
  hi  <- 1 - eps
  
  # candidate set: roots that lie in (pi_min, 1)
  cand <- c(r1, r2)
  cand <- cand[is.finite(cand) & (cand > lo) & (cand < hi)]
  
  # if no root is admissible (rare), fall back to interior projection of A/n
  if (length(cand) == 0) {
    p0 <- A / n
    return(min(max(p0, lo), hi))
  }
  
  # evaluate objective and pick best root
  f <- function(pi) A*log(pi) + B*log(1 - pi) + bw*log(pi - pi_min)
  vals <- sapply(cand, f)
  cand[which.max(vals)]
}


# --- ZIP M-step with optional barrier on pi ------------------------------

zip_mstep_barrier <- function(x, gamma, use_barrier = FALSE, bw = 0, pi_min = 0.005) {
  n <- length(x)
  A <- sum(1 - gamma)
  B <- sum(gamma)
  
  # lambda update (same as usual)
  lambda_new <- sum(gamma * x) / sum(gamma)
  
  # pi update
  if (!use_barrier) {
    pi_new <- A / n
  } else {
    pi_new <- zip_mstep_pi_barrier(A, B, bw = bw, pi_min = pi_min)
  }
  
  list(pi = pi_new, lambda = lambda_new)
}

# ZIP barrier-EM -----------------------------------------------------------
# - r is fixed (default r=1 => standard posterior in E-step)
# - barrier schedule bw_init -> bw_end over n_steps
# - uses your existing: zip_estep_annealed(), zip_mstep_barrier()

zip_barrier <- function(x,
                        theta_init,
                        bw_init = 5,
                        bw_end  = 0.05,
                        pi_min  = 0.005,
                        n_steps = 50,
                        tol = 1e-10,
                        max_iter_inner = 200) {
  
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))
  
  theta <- theta_init
  
  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (s in seq_len(n_steps)) {
    bw <- bw_grid[s]
    
    for (it in seq_len(max_iter_inner)) {
      pi_old     <- theta$pi
      lambda_old <- theta$lambda
      
      # E-step (r=1 고정)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = 1)
      
      # M-step (barrier on pi)
      theta <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )
      
      d_pi  <- abs(theta$pi - pi_old)
      d_lam <- abs(theta$lambda - lambda_old)
      
      trace <- rbind(trace, data.frame(
        step = s, inner = it, bw = bw,
        pi = theta$pi, lambda = theta$lambda,
        d_pi = d_pi, d_lam = d_lam
      ))
      
      if (max(d_pi, d_lam) < tol) break
    }
  }
  
  list(theta = theta, trace = trace, bw_grid = bw_grid)
}


# --- you provided (keep as-is) -------------------------------------------
zip_dkl_std_post <- function(x, theta0, theta1) {
  z0 <- (x == 0L)
  if (!any(z0)) return(0)
  
  g0 <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1 <- zip_estep_annealed(x[z0], pi = theta1$pi, lambda = theta1$lambda, r = 1)
  
  sum(g0 * (log(g0) - log(g1)) + (1 - g0) * (log(1 - g0) - log(1 - g1)))
}


zip_delta_dkl <- function(x, theta0, theta1, r) {
  z0 <- (x == 0L)
  if (!any(z0)) return(0)
  
  g_r0 <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = r)
  g0   <- zip_estep_annealed(x[z0], pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1   <- zip_estep_annealed(x[z0], pi = theta1$pi, lambda = theta1$lambda, r = 1)
  
  sum(g_r0 * (log(g0) - log(g1)) + (1 - g_r0) * (log(1 - g0) - log(1 - g1)))
}

zip_delta_dkl_debug <- function(x, theta0, theta1, r, eps = 0, verbose = TRUE) {
  z0 <- (x == 0L)
  if (!any(z0)) return(list(value = 0, info = "no zeros"))

  x0 <- x[z0]

  g_r0 <- zip_estep_annealed(x0, pi = theta0$pi, lambda = theta0$lambda, r = r)
  g0   <- zip_estep_annealed(x0, pi = theta0$pi, lambda = theta0$lambda, r = 1)
  g1   <- zip_estep_annealed(x0, pi = theta1$pi, lambda = theta1$lambda, r = 1)

  # (선택) eps>0이면 클리핑해서 NaN을 없애고 값은 계속 계산 가능
  if (eps > 0) {
    clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)
    g_r0 <- clip01(g_r0); g0 <- clip01(g0); g1 <- clip01(g1)
  }

  # 로그 계산 전 점검
  bad <- list(
    g_r0_nonfinite = which(!is.finite(g_r0)),
    g0_nonfinite   = which(!is.finite(g0)),
    g1_nonfinite   = which(!is.finite(g1)),
    g0_le0         = which(g0 <= 0),
    g0_ge1         = which(g0 >= 1),
    g1_le0         = which(g1 <= 0),
    g1_ge1         = which(g1 >= 1),
    one_minus_g0_le0 = which((1 - g0) <= 0),
    one_minus_g1_le0 = which((1 - g1) <= 0)
  )

  # 실제로 log가 Inf/NaN을 만드는지 확인
  lg0  <- log(g0)
  lg1  <- log(g1)
  l1g0 <- log(1 - g0)
  l1g1 <- log(1 - g1)

  bad$log_g0_nonfinite   <- which(!is.finite(lg0))
  bad$log_g1_nonfinite   <- which(!is.finite(lg1))
  bad$log1mg0_nonfinite  <- which(!is.finite(l1g0))
  bad$log1mg1_nonfinite  <- which(!is.finite(l1g1))

  # 항별로도 NaN 생기는지 체크
  term1 <- g_r0 * (lg0 - lg1)
  term2 <- (1 - g_r0) * (l1g0 - l1g1)

  bad$term1_nonfinite <- which(!is.finite(term1))
  bad$term2_nonfinite <- which(!is.finite(term2))

  value <- sum(term1 + term2)

  if (verbose) {
    cat("---- zip_delta_dkl_debug ----\n")
    cat("n0 =", length(x0), "\n")
    cat("ranges:\n")
    cat("  g_r0:", sprintf("[%.3g, %.3g]", min(g_r0), max(g_r0)), "\n")
    cat("  g0  :", sprintf("[%.3g, %.3g]", min(g0), max(g0)), "\n")
    cat("  g1  :", sprintf("[%.3g, %.3g]", min(g1), max(g1)), "\n")
    cat("value =", value, "\n")

    # 문제가 있는 인덱스만 요약 출력
    show_idx <- unique(c(
      bad$g_r0_nonfinite, bad$g0_nonfinite, bad$g1_nonfinite,
      bad$g0_le0, bad$g0_ge1, bad$g1_le0, bad$g1_ge1,
      bad$log_g0_nonfinite, bad$log_g1_nonfinite,
      bad$log1mg0_nonfinite, bad$log1mg1_nonfinite,
      bad$term1_nonfinite, bad$term2_nonfinite
    ))
    show_idx <- show_idx[!is.na(show_idx)]

    if (length(show_idx) == 0) {
      cat("No non-finite issues detected.\n")
    } else {
      show_idx <- head(show_idx, 10)  # 너무 길면 10개만
      cat("Problem indices (first up to 10):", paste(show_idx, collapse=", "), "\n")
      df <- data.frame(
        i = show_idx,
        x0 = x0[show_idx],
        g_r0 = g_r0[show_idx],
        g0   = g0[show_idx],
        g1   = g1[show_idx],
        logg0 = lg0[show_idx],
        logg1 = lg1[show_idx],
        log1mg0 = l1g0[show_idx],
        log1mg1 = l1g1[show_idx],
        term1 = term1[show_idx],
        term2 = term2[show_idx]
      )
      print(df)
      cat("Bad sets sizes:\n")
      print(sapply(bad, length))
    }
    cat("-----------------------------\n")
  }

  list(value = value, bad = bad)
}


# DHEM with fixed schedules ------------------------------------------------
# - r_t : exp(seq(log(r_init), 0, length.out = n_steps))  (monotone up to 1)
# - bw_t: exp(seq(log(bw_init), log(bw_end), length.out = n_steps)) (monotone down if bw_end<bw_init)
# - At each step: E-step with annealed posterior (r_t),
#                 M-step with barrier on pi using bw_t

zip_dhem <- function(x,
                     theta_init,
                     r_init = 0.2,
                     bw_init = 5,
                     bw_end  = 0.05,
                     pi_min  = 0.005,
                     n_steps = 50,
                     tol = 1e-10,
                     max_iter_inner = 200) {
  
  # fixed schedules
  r_grid  <- exp(seq(log(r_init), 0, length.out = n_steps))                 # -> 1
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = n_steps))      # bw_init -> bw_end
  
  theta <- theta_init
  
  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    r      = numeric(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  
  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]
    bw <- bw_grid[s]
    
    for (it in seq_len(max_iter_inner)) {
      pi_old     <- theta$pi
      lambda_old <- theta$lambda
      
      # E-step (annealed with r)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = r)
      
      # M-step (barrier on pi with bw)
      theta <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )
      
      d_pi  <- abs(theta$pi - pi_old)
      d_lam <- abs(theta$lambda - lambda_old)
      
      trace <- rbind(trace, data.frame(
        step = s, inner = it, r = r, bw = bw,
        pi = theta$pi, lambda = theta$lambda,
        d_pi = d_pi, d_lam = d_lam,
        stringsAsFactors = FALSE
      ))
      
      if (max(d_pi, d_lam) < tol) break
    }
  }
  
  list(theta = theta, trace = trace, r_grid = r_grid, bw_grid = bw_grid)
}


zip_adaptive_dhem <- function(x,
                              theta_init,
                              r_init  = 0.2,
                              bw_init = 1e-6,
                              pi_min  = 5e-6,
                              n_steps = 50,
                              eta = 0.1,
                              bw_rate = 0.1,
                              tol = 1e-10,
                              max_inner = 500,
                              verbose = TRUE
                            ) {
  
  r_grid   <- exp(seq(log(r_init), 0, length.out = n_steps))
  bw = bw_init
  theta <- theta_init
  
  trace <- data.frame(
    step   = integer(0),
    inner  = integer(0),
    r      = numeric(0),
    bw     = numeric(0),
    pi     = numeric(0),
    lambda = numeric(0),
    d_pi   = numeric(0),
    d_lam  = numeric(0),
    stringsAsFactors = FALSE
  )
  theta0 = list(pi=theta_init$pi,lambda=theta_init$lambda)

  for (s in seq_len(n_steps)) {
    r  <- r_grid[s]

    theta_temp = theta0
    for (it in seq_len(max_inner)) {
      pi_old     <- theta_temp$pi
      lambda_old <- theta_temp$lambda
      # E-step (annealed with r)
      gamma <- zip_estep_annealed(x, pi = pi_old, lambda = lambda_old, r = r)
      # M-step (barrier on pi with bw)
      theta1 <- zip_mstep_barrier(
        x, gamma,
        use_barrier = TRUE,
        bw = bw,
        pi_min = pi_min
      )
      # print(theta1$pi)

      d_pi  <- abs(theta1$pi - theta_temp$pi )
      d_lam <- abs(theta1$lambda - theta_temp$lambda)
      
      if (max(d_pi, d_lam) < tol) {break
        }else{
        theta_temp=theta1
      }
    }

    dDKL   <- zip_delta_dkl(x, theta0 = theta0, theta1 = theta1, r = r)
    DKLstd <- zip_dkl_std_post(x, theta0 = theta0, theta1 = theta1)
    dB <- log(theta1$pi - pi_min) - log(theta0$pi - pi_min)

    Accept = FALSE
    if(is.na(dDKL-bw*dB)||dDKL-bw*dB<0){
      # Acceptance test start
      # Acc1
      if(dDKL<=eta*DKLstd){
        next
      }
      # Acc2
      if(eta*DKLstd<bw*abs(dB)){
        bw = min(bw,eta*DKLstd/abs(dB))
        Accept = FALSE
      }
    }else{
      Accept=TRUE
    }

    # parameter update 
    if(Accept) theta0 = theta1
    #print(theta0$pi)
    trace <- rbind(trace, data.frame(
        step = s, inner = it, r = r, bw = bw,
        pi = theta0$pi, lambda = theta0$lambda,
        d_pi = d_pi, d_lam = d_lam,
        stringsAsFactors = FALSE
      ))
      
  }
  
  list(theta = theta0, bw = bw)
}







# data generate -----------------------------------------------------------------

n <- 1e+4
pi_true  <- 0.99
lambda_true <- 0.3
x <- rzip(n, pi = pi_true, lambda = lambda_true)

pi_init_value = 0.7
lambda_init_value = 1
theta_init_value = list(pi = pi_init_value, lambda = lambda_init_value)
# lambda_init = compute_lambda_init(x,pi_init_value)
r_init_value = 0.1
p_min_value = 0.5
n_step_value = 100
bw_init_value = 1e+1
bw_end_value = 1e-5

res_em <- zip_em(x, pi_init = pi_init_value)
res_adapDhem <- zip_adaptive_dhem(
  x = x,
  theta_init = theta_init_value,
  r_init  = r_init_value,
  bw_init = bw_init_value,
  pi_min  = p_min_value,
  n_steps = n_step_value,
  bw_rate = 0.9,
  tol = 1e-10,
  max_inner = 1e+4,
  verbose = FALSE
)
res_em$theta
res_adapDhem$theta


res_daem <- zip_daem(x,theta_init = theta_init_value,r_init = r_init_value,n_steps  = n_step_value)
res_barrier <- zip_barrier( x,theta_init = theta_init_value,n_steps = n_step_value,
  bw_init = bw_init_value,bw_end  = bw_end_value,pi_min  = p_min_value)
res_dhem <- zip_dhem( x,theta_init = theta_init_value,r_init = r_init_value, bw_init = bw_init_value,
    bw_end  = bw_end_value,  pi_min = p_min_value,  n_steps = n_step_value)
res_adapDhem <- zip_adaptive_dhem(
  x = x,
  theta_init = theta_init_value,
  r_init  = r_init_value,
  bw_init = bw_init_value,
  pi_min  = p_min_value,
  n_steps = n_step_value,
  bw_rate = 0.9,
  tol = 1e-10,
  max_inner = 1e+4,
  verbose = FALSE
)

res_em$theta
# res_daem$theta
# res_adapDaem$theta
res_barrier$theta
res_dhem$theta
res_adapDhem$theta



sum(x==0)

tail(res_adapDhem)

pi_err_em        <- log_err(res_em$theta$pi,        pi_true)
pi_err_daem      <- log_err(res_daem$theta$pi,      pi_true)
# pi_err_adapDaem  <- log_err(res_adapDaem$theta$pi,  pi_true)
pi_err_barrier   <- log_err(res_barrier$theta$pi,   pi_true)
pi_err_dhem      <- log_err(res_dhem$theta$pi,   pi_true)
pi_err_adapdhem  <- log_err(res_adapDhem$theta$pi,   pi_true)

lambda_err_em        <- abs_err(res_em$theta$lambda,        lambda_true)
lambda_err_daem      <- abs_err(res_daem$theta$lambda,      lambda_true)
# lambda_err_adapDaem  <- abs_err(res_adapDaem$theta$lambda,  lambda_true)
lambda_err_barrier   <- abs_err(res_barrier$theta$lambda,   lambda_true)
lambda_err_dhem      <- abs_err(res_dhem$theta$lambda   ,   lambda_true)
lambda_err_adapdhem  <- abs_err(res_adapDhem$theta$lambda,  lambda_true)

c(
  EM           = pi_err_em,
  DAEM         = pi_err_daem,
  # AdaptiveDAEM = pi_err_adapDaem,
  BarrierEM    = pi_err_barrier,
  DHEM = pi_err_dhem,
  AdaptiveDHEM = pi_err_adapdhem
)


## Simulation 

# ----------------------------
# Simulation
# ----------------------------
M <- 100
sim_df <- data.frame(
  simulationNum = integer(0),
  method        = character(0),
  pi_err        = numeric(0),
  lambda_err    = numeric(0),
  stringsAsFactors = FALSE
)

M=100
t0 <- Sys.time()
sim_df <- data.frame()

for (i in 1:M) {
  
  if (i %% 10 == 0 || i == 1) {
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- dt / i
    eta  <- rate * (M - i)
    message(sprintf(
      "Simulation %d/%d (%.1f%%) | elapsed=%.1fs | ETA=%.1fs",
      i, M, 100*i/M, dt, eta
    ))
  }
  
n <- 1e+4
pi_true  <- 0.99
lambda_true <- 0.3
x <- rzip(n, pi = pi_true, lambda = lambda_true)

pi_init_value = 0.7
lambda_init_value = 1
theta_init_value = list(pi = pi_init_value, lambda = lambda_init_value)
# lambda_init = compute_lambda_init(x,pi_init_value)
r_init_value = 0.1
p_min_value = 0.5
n_step_value = 200
bw_init_value = 1e-1
bw_end_value = 1e-5
  
  ## helper: 실패해도 loop 계속 + succ 생성
  safe_fit <- function(expr) {
    out <- tryCatch(expr, error = function(e) NULL)
    if (is.null(out) || is.null(out$theta) || is.null(out$theta$pi) || is.null(out$theta$lambda)) {
      return(list(succ = 0L, theta = list(pi = NA_real_, lambda = NA_real_)))
    } else {
      return(list(succ = 1L, theta = out$theta))
    }
  }
  
  ## fit (각 method별로 succ/theta 확보)
  res_em       <- safe_fit(zip_em(x, pi_init = pi_init_value))
  res_daem     <- safe_fit(zip_daem(x, theta_init = theta_init_value, r_init = r_init_value, n_steps = n_step_value))
  # res_adapDaem <- safe_fit(zip_adaptive_daem(x, theta_init = theta_init_value, r_init = r_init_value,
  #                                            n_steps = n_step_value, r_switch = 0, eta = 0.5))
  res_barrier  <- safe_fit(zip_barrier(x, theta_init = theta_init_value, n_steps = n_step_value,
                                       bw_init = bw_init_value, bw_end = bw_end_value, pi_min = p_min_value))
  res_dhem     <- safe_fit(zip_dhem(x, theta_init = theta_init_value, r_init = r_init_value,
                                    bw_init = bw_init_value, bw_end = bw_end_value, pi_min = p_min_value,
                                    n_steps = n_step_value))
  res_adapDhem <- safe_fit(zip_adaptive_dhem(
  x = x,
  theta_init = theta_init_value,
  r_init  = r_init_value,
  bw_init = bw_init_value,
  pi_min  = p_min_value,
  n_steps = n_step_value,
  bw_rate = 0.9,
  tol = 1e-10,
  max_inner = 1e+4,
  verbose = FALSE
  ))
  
  ## collect (long format)
  sim_df <- rbind(
    sim_df,
    data.frame(
      simulationNum = i,
      method = "EM",
      succ = res_em$succ,
      pi_err = ifelse(res_em$succ == 1L, log_err(res_em$theta$pi, pi_true), NA_real_),
      lambda_err = ifelse(res_em$succ == 1L, log_err(res_em$theta$lambda, lambda_true), NA_real_),
      stringsAsFactors = FALSE
    ),
    data.frame(
      simulationNum = i,
      method = "DAEM",
      succ = res_daem$succ,
      pi_err = ifelse(res_daem$succ == 1L, err(res_daem$theta$pi, pi_true), NA_real_),
      lambda_err = ifelse(res_daem$succ == 1L, err(res_daem$theta$lambda, lambda_true), NA_real_),
      stringsAsFactors = FALSE
    ),
    data.frame(
      simulationNum = i,
      method = "Barrier",
      succ = res_barrier$succ,
      pi_err = ifelse(res_barrier$succ == 1L, err(res_barrier$theta$pi, pi_true), NA_real_),
      lambda_err = ifelse(res_barrier$succ == 1L, err(res_barrier$theta$lambda, lambda_true), NA_real_),
      stringsAsFactors = FALSE
    ),
    data.frame(
      simulationNum = i,
      method = "DHEM",
      succ = res_dhem$succ,
      pi_err = ifelse(res_dhem$succ == 1L, err(res_dhem$theta$pi, pi_true), NA_real_),
      lambda_err = ifelse(res_dhem$succ == 1L, err(res_dhem$theta$lambda, lambda_true), NA_real_),
      stringsAsFactors = FALSE
    ),
    data.frame(
      simulationNum = i,
      method = "AdaptiveDHEM",
      succ = res_adapDhem$succ,
      pi_err = ifelse(res_adapDhem$succ == 1L, err(res_adapDhem$theta$pi, pi_true), NA_real_),
      lambda_err = ifelse(res_adapDhem$succ == 1L, err(res_adapDhem$theta$lambda, lambda_true), NA_real_),
      stringsAsFactors = FALSE
    )
  )
}

########## Table 

library(dplyr)
library(tidyr)
library(ggplot2)

# 1) method별 err mean/sd 테이블
library(dplyr)

err_tbl_disp <- sim_df %>%
  group_by(method) %>%
  summarise(
    pi_mean     = mean(pi_err, na.rm = TRUE),
    pi_sigma    = sd(pi_err,   na.rm = TRUE),
    lambda_mean = mean(lambda_err, na.rm = TRUE),
    lambda_sigma= sd(lambda_err,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pi      = sprintf("%.3f (%.3f)", pi_mean, pi_sigma),
    lambda  = sprintf("%.3f (%.3f)", lambda_mean, lambda_sigma)
  ) %>%
  select(method, pi, lambda)

err_tbl_disp

# 2) boxplot + facet(pi/lambda)

plot_df <- sim_df %>%
  pivot_longer(cols = c(pi_err, lambda_err),
               names_to = "param", values_to = "err") %>%
  mutate(param = recode(param, pi_err = "pi", lambda_err = "lambda"))

ggplot(plot_df, aes(x = method, y = err)) +
  geom_boxplot(outlier.alpha = 0.5) +
  facet_wrap(~ param, scales = "free_y",
             labeller = labeller(param = label_parsed)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=16),
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 20)  ,
    axis.title.x = element_blank()   # <- "method"만 제거
  )














###############

library(dplyr)
library(tidyr)

# 1) 모든 방법이 성공한 simulation만 필터
all_ok_ids <- sim_df %>%
  group_by(simulationNum) %>%
  summarise(all_ok = all(succ == 1L), .groups = "drop") %>%
  filter(all_ok) %>%
  pull(simulationNum)

df_ok <- sim_df %>%
  filter(simulationNum %in% all_ok_ids)

# 2) EM error를 simulationNum별로 붙이고, EM-method 차이(delta) 계산
delta_df <- df_ok %>%
  select(simulationNum, method, pi_err, lambda_err) %>%
  pivot_longer(cols = c(pi_err, lambda_err),
               names_to = "para", values_to = "err") %>%
  group_by(simulationNum, para) %>%
  mutate(err_em = err[method == "EM"]) %>%
  ungroup() %>%
  filter(method != "EM") %>%                # x축에 EM 제외 (원하면 유지 가능)
  mutate(delta = err - err_em  ) %>%          # 질문의 "EM - method"
  select(simulationNum, method, para, delta)

# para 라벨 정리 (facet용)
delta_df <- delta_df %>%
  mutate(para = recode(para,
                       pi_err = "pi",
                       lambda_err = "lambda"))


library(ggplot2)

ggplot(delta_df, aes(x = method, y = delta)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(
    ~ para, scales = "free_y",
    labeller = as_labeller(c(pi = "π", lambda = "λ"))
  ) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(
    title = "Difference error (EM - method)",
    subtitle = "Difference < 0 means smaller error than EM",
    x = "",
    y = "Difference"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
strip.text = element_text(size = 16))


##########



M <- 200  # 원하는 반복 횟수

em_list <- vector("list", M)

for (m in seq_len(M)) {
  x <- rzip(n, pi = pi_true, lambda = lambda_true)
  
  # zip_em의 인자/반환 구조에 맞게 수정하세요
  res_em <- zip_em(x, pi_init = pi_init_value, lambda_init = lambda_init_value)
  
  # (A) 가장 흔한 케이스: res_em$pi, res_em$lambda
  # (B) 혹은 res_em$theta$pi, res_em$theta$lambda
  pi_hat <- if (!is.null(res_em$pi)) res_em$pi else res_em$theta$pi
  lam_hat <- if (!is.null(res_em$lambda)) res_em$lambda else res_em$theta$lambda
  
  em_list[[m]] <- data.frame(sim = m, pi_hat = pi_hat, lambda_hat = lam_hat)
}

em_df <- bind_rows(em_list)


p_pi <- ggplot(em_df, aes(x = 1, y = pi_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = pi_true, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "π (EM estimate)", x = NULL, y = "") +
    coord_cartesian(ylim = c(0.85,1))+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=15),
      )

p_lam <- ggplot(em_df, aes(x = 1, y = lambda_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = lambda_true, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "λ (EM estimate)", x = NULL, y = "") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=15))

(p_pi | p_lam) +
  plot_annotation(
    title = sprintf("ZIP- standard EM estimates over %d simulations", M),
        subtitle = sprintf("n=%d, true π=%.3f, true λ=%.3f ",
                       n, pi_true, lambda_true))


###############

M <- 200
adap_list <- vector("list", M)

for (m in seq_len(M)) {
  x <- rzip(n, pi = pi_true, lambda = lambda_true)

  res <- zip_adaptive_dhem(
    x = x,
    theta_init = theta_init_value,
    r_init  = r_init_value,
    bw_init = bw_init_value,
    pi_min  = p_min_value,
    n_steps = n_step_value,
    bw_rate = 0.9,
    tol = 1e-10,
    max_inner = 1e+4,
    verbose = FALSE
  )

  # 최종 추정치 추출 (반환 구조에 맞춰 조정)
  pi_hat  <- res$theta$pi
  lam_hat <- res$theta$lambda

  # 실패(NA/Inf) 방어
  adap_list[[m]] <- data.frame(
    sim = m,
    pi_hat = pi_hat,
    lambda_hat = lam_hat,
    bw_final = if (!is.null(res$bw)) res$bw else NA_real_
  )
}

adap_df <- bind_rows(adap_list) %>%
  filter(is.finite(pi_hat), is.finite(lambda_hat))

# ---- plots (patchwork) ----
p_pi <- ggplot(adap_df, aes(x = 1, y = pi_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = pi_true, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "π (Adaptive DHEM estimate)", x = NULL, y = "") +
    coord_cartesian(ylim = c(0.85,1))+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

p_lam <- ggplot(adap_df, aes(x = 1, y = lambda_hat)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_hline(yintercept = lambda_true, color = "blue", linetype = 2, linewidth = 0.8) +
  labs(title = "λ (Adaptive DHEM estimate)", x = NULL, y = "") +
  coord_cartesian(ylim=c(0,0.8))+
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 15))

(p_pi | p_lam) +
  plot_annotation(
    title = sprintf("ZIP - Adaptive DHEM estimates over %d simulations", M),
    subtitle = sprintf("n=%d, true π=%.3f, true λ=%.3f ",
                       n, pi_true, lambda_true)
  )



#######



