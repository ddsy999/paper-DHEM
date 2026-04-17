source("zip_functions.R")





# data generate -----------------------------------------------------------------

n <- 1e+4
pi_true  <- 0.99
lambda_true <- 0.3
x <- rzip(n, pi = pi_true, lambda = lambda_true)

pi_init_value = 0.7
lambda_init_value = 1
theta_init_value = list(pi = pi_init_value, lambda = lambda_init_value)
r_init_value = 0.1
p_min_value = 0.5
n_step_value = 100
bw_end_value = 1e-5

bw_init_value = zip_bw_init_pi(x,pi_init_value,lambda_init_value,r_init_value,p_min_value,tau=0.1)

## Simulation 

# ----------------------------
# Simulation
# ----------------------------
M <- 100 #200
sim_df <- data.frame(
  simulationNum = integer(0),
  method        = character(0),
  pi_err        = numeric(0),
  lambda_err    = numeric(0),
  stringsAsFactors = FALSE
)

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
lambda_true <- 0.4
x <- rzip(n, pi = pi_true, lambda = lambda_true)

pi_init_value = 0.7
lambda_init_value = 1
theta_init_value = list(pi = pi_init_value, lambda = lambda_init_value)
# lambda_init = compute_lambda_init(x,pi_init_value)
r_init_value = 0.1
p_min_value = 0.5
n_step_value = 100
#bw_init_value = 1e-1 
bw_init_value = min(zip_bw_init_pi(x,pi_init_value,lambda_init_value,r_init_value,p_min_value,tau=0.1))

bw_end_value = 1e-8
  
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



M <- 100  # 원하는 반복 횟수

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
        subtitle = sprintf("n=%d, true π=%.3f, true λ=%.4f ",
                       n, pi_true, lambda_true))


###############

M <- 100
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
    subtitle = sprintf("n=%d, true π=%.3f, true λ=%.4f ",
                       n, pi_true, lambda_true)
  )



#######



