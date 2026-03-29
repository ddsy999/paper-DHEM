
source("gmm_functions.R")

## ------------------------------------------------------------
## Simulation
## ------------------------------------------------------------

theta_true <- list(
  pi = c(0.2,0.3,0.5),
  mu = list(c(1,1,2), c(5,3,0.5), c(2,6,5)),
  Sigma = list(
    matrix(c(1,0,0,0,1,0,0,0,1),3,3),
    matrix(c(2,0,0,0,0.1,0,0,0,0.5),3,3),
    matrix(c(0.5,0,0,0,0.1,0,0,0,2),3,3)
  )
)

theta_true <- list(
  pi = c(0.2,0.3,0.5),
  mu = list(c(1,1,2), c(5,3,0.5), c(2,6,5)),
  Sigma = list(
    # diag = (1, 1, 1)
    matrix(c(
      1,    0.30, -0.20,
      0.30, 1,     0.10,
      -0.20, 0.10,  1
    ), 3, 3, byrow = TRUE),
    
    # diag = (2, 0.1, 0.5)
    matrix(c(
      2,     0.18, -0.25,
      0.18,  0.10,  0.06,
      -0.25,  0.06,  0.50
    ), 3, 3, byrow = TRUE),
    
    # diag = (0.5, 0.1, 2)
    matrix(c(
      0.50, -0.08,  0.22,
      -0.08,  0.10, -0.12,
      0.22, -0.12,  2
    ), 3, 3, byrow = TRUE)
  )
)


n=100
K=3
theta_true = theta_true
tol_limit = 1e-6
r_init = 0.1
n_steps = 200
bw_init = 1e-1
bw_end = 1e-3
delta =  10 #10 #200
max_iter= 1000
eta = 0.02

## ------------------------------------------------------------
## Simulation M 
## ------------------------------------------------------------


M=100
sim_list <- lapply(seq_len(M), function(m) {
  if (m == 1 || m %% 10 == 0) message(sprintf("Simulation %d/%d (%.1f%%)", m, M, 100*m/M))
  one_sim_summary(simulation_num = m, n = n, K = K, theta_true = theta_true,tol_limit=tol_limit,
                  r_init=r_init,n_steps=n_steps,bw_init=bw_init,bw_end=bw_end,delta=delta,
                  max_iter=max_iter,eta=eta)
})
sim_all <- do.call(rbind, sim_list)


#########

library(dplyr)

err_summary <- sim_all %>%
  filter(!is.na(err)) %>%        # 실패(run) 제거
  group_by(method, para) %>%
  summarise(
    mean_err = mean(err),
    sd_err   = sd(err),
    n        = n(),
    .groups = "drop"
  )

err_table <- sim_all %>%
  filter(succ == 1) %>%
  group_by(method, para) %>%
  summarise(
    value = sprintf("%.3f (%.3f)", mean(err), sd(err)),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = para,
    values_from = value
  )

err_table



#########

ggplot(sim_all %>% filter(!is.na(err)),
       aes(method, err)) +
  geom_boxplot() +
  facet_wrap(
    ~ para,
    scales = "free_y",
    labeller = as_labeller(
      c(mu = "mu", Sigma = "Sigma", pi = "pi"),
      label_parsed
    )
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=16),
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank()
  )





sim_all[which(sim_all[,"succ"]==0),]


success_table <- sim_all %>%
  group_by(method) %>%
  summarise(
    n_total = n(),
    n_success = sum(succ == 1, na.rm = TRUE),
    success_rate = n_success / n_total,
    .groups = "drop"
  ) %>%
  arrange(desc(success_rate))

success_table

