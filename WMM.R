

source("wmm_functions.R")

#### ---------------------------------------------------------------------------
# Data Load
#### ---------------------------------------------------------------------------
getwd()
df <- read.table("DATA\\Aarest_data.txt", header = TRUE)


#### ---------------------------------------------------------------------------
# Hyper-parameter (set here first)
#### ---------------------------------------------------------------------------
maxGEMiter = 1e+6
nsteps     <- 100
r_init     <- 0.1
r_end      <- 1
bw_init    <- 1e-1
bw_end     <- 1e-8   # 필요시 조정
eta        <- 0.1           # adaptive DHEM 전용
errtol     = 1e-10


#### ---------------------------------------------------------------------------
# Init-parameter
#### ---------------------------------------------------------------------------
K <- 3
pi_init   <- rep(1 / K, K)
beta_init <- c(0.5, 1, 2)   
lambda_init <- wmm_lambda_init(df$time,df$event, beta_init,ratio1=0.3,ratio3=0.7)
theta_init = list(beta=beta_init,pi=pi_init,lambda=lambda_init)

#### ---------------------------------------------------------------------------
# Train (one run each)
#### ---------------------------------------------------------------------------
verbose=FALSE
# 1) EM (standard EM: r=1, bw=0) 
fit_EM <- wmm_EM(df,theta=theta_init,maxGEMiter = maxGEMiter,tol=errtol,verbose=verbose)
# 2) DAEM (bw=0, r schedule)
fit_DAEM <- wmm_DAEM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps = nsteps,r_init = r_init,r_end=r_end,tol=errtol,verbose=verbose)
# 3) Barrier method (r=1, bw schedule)
fit_BM <- wmm_BM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
# 4) DHEM (r schedule + bw schedule)
fit_DHEM <- wmm_DHEM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
# 5) Adaptive DHEM (r schedule + adaptive bw control)
fit_adapDHEM <- wmm_DHEM_adaptive(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,eta=eta,tol = errtol,verbose=verbose)



df_EM       <- make_result_df(fit_EM, "EM")
df_DAEM     <- make_result_df(fit_DAEM, "DAEM")
df_BM       <- make_result_df(fit_BM, "BM")
df_DHEM     <- make_result_df(fit_DHEM, "DHEM")
df_adapDHEM <- make_result_df(fit_adapDHEM, "AdaptiveDHEM")



#### ---------------------------------------------------------------------------
# Table
#### ---------------------------------------------------------------------------
result_table <- do.call(rbind, lapply(
  list(df_EM, df_DAEM, df_BM, df_DHEM, df_adapDHEM),
  best_row
))

result_table

tail(df_DHEM,30)
#### ---------------------------------------------------------------------------
# Figure
#### ---------------------------------------------------------------------------
(plot_beta_trace(df_DAEM, title = expression(DAEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_DAEM, title = expression(DAEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "DAEM") &  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace_barrier(df_BM, title = expression(Barrier*" "*method * ": " * beta * " trace"))|plot_dQbeta_trace_barrier(df_BM, title = expression(DHEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "Barrier method")&  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace(df_DHEM, title = expression(DHEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_DHEM, title = expression(DHEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "DHEM")&  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace(df_adapDHEM, title = expression(adapDHEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_adapDHEM, title = expression(adapDHEM* ": " * nabla~beta * " trace")))+
  plot_annotation(title= "Adaptive DHEM")&  theme(plot.title = element_text(size = 20, hjust = 0.5))


tail(df_adapDHEM,1)
tail(df_BM)



 do.call(rbind, lapply(
  list(df_EM, df_DAEM, df_BM, df_DHEM, df_adapDHEM),
  function(x) tail(x,1)
))


