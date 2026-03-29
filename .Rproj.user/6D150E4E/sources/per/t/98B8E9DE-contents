
diffB_onlyB = function(beta,event_vec,time_vec,latentZ_mat,j){
  sum(latentZ_mat[,j]*event_vec)/beta + 
    sum(latentZ_mat[,j]*event_vec*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)*sum(latentZ_mat[,j]*(time_vec^beta)*log(time_vec))/sum(latentZ_mat[,j]*(time_vec^beta))
}

barrierFunc_1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1)+(1/beta -1/(1-beta))*(bw)
  return(result)
}

barrierFunc_3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3)+bw*(1/(beta-1))
  return(result)
}

barrier_beta1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1))

  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps = 1e-12
  if(bw==0){
    while(!isna_diffbeta(maxRange)){
      maxRange=maxRange + min(maxRange*1.01,10)
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=1)*
      diffB_onlyB(1e-12,event_vec,time_vec,latentZ_mat, j=1)<0) break
    }
      result <- uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1),
      interval = c(0,maxRange),tol=1e-10)
      return(result$root)
  }

  result <- uniroot(function(beta) barrierFunc_1(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(eps,1-eps),tol=1e-10)
  return(result$root)
}

barrier_beta3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3))
  
  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps= 1e-3
  if(bw==0){
    while(!isna_diffbeta(maxRange)&&!isna_diffbeta(eps)){
      maxRange=maxRange*1.01
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
      diffB_onlyB(eps,event_vec,time_vec,latentZ_mat, j=3)<0) break
    }
    result = uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3),
    interval = c(eps, maxRange),tol=1e-10)
    return(result$root)
  }

  while(!isna_diffbeta(maxRange)){
    maxRange=maxRange + min(maxRange*1.01,10)
    if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
    diffB_onlyB(1,event_vec,time_vec,latentZ_mat, j=3)<0) break
  }
  result = uniroot(function(beta) barrierFunc_3(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}

barrier_safe_wrapper1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta1(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

barrier_safe_wrapper3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta3(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

weibull_estep_annealed <- function(df, pi, lambda, beta, r = 1) {
  t     <- as.numeric(df$time)
  event <- as.numeric(df$event)
  
  n <- length(t)
  K <- length(pi)
  
  logt <- log(t)
  log_resp <- matrix(NA_real_, nrow = n, ncol = K)
  
  for (k in 1:K) {
    # Weibull component k:
    # f_k(t) = lambda_k * beta_k * t^(beta_k-1) * exp(-lambda_k * t^beta_k)
    # S_k(t) = exp(-lambda_k * t^beta_k)
    #
    # log L_ik = event_i * log f_k(t_i) + (1-event_i) * log S_k(t_i)
    #         = event_i*(log lambda_k + log beta_k + (beta_k-1)log t_i) - lambda_k t_i^beta_k
    logLik_ik <- event * (log(lambda[k]) + log(beta[k]) + (beta[k] - 1) * logt) -
      lambda[k] * (t ^ beta[k])
    
    # annealed responsibilities: gamma_ik ∝ exp( r * (log pi_k + logLik_ik) )
    log_resp[, k] <- r * (log(pi[k]) + logLik_ik)
  }
  
  # numeric stabilization and normalization across k
  rowmax <- apply(log_resp, 1, max)
  w <- exp(log_resp - rowmax)
  gamma <- w / rowSums(w)
  
  gamma
}

wmm_delta_DKL <- function(df, theta0, theta1, r) {

  # standard posterior (r=1) for theta0, theta1
  gamma0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  gamma1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)

  # annealed posterior weight (r=r) under theta0
  w <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = r)

  # per-i contribution: sum_k w_ik * (log gamma0_ik - log gamma1_ik)
  contrib <- rowSums(w * (log(gamma0) - log(gamma1)))

  mean(contrib)
}

wmm_DKL <- function(df, theta0, theta1) {
  g0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  g1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)
  mean(rowSums(g0 * (log(g0) - log(g1))))
}

wmm_bar_value <- function(beta) {
  log(beta[1]) + log(1 - beta[1]) + log(beta[3] - 1)
}

wmm_bar_diff <- function(theta_new, theta_old) {
  wmm_bar_value(theta_new$beta) - wmm_bar_value(theta_old$beta)
}

wmm_lambda_init <-function(time_vec,event_vec,beta_vec,ratio1,ratio3){
  surv_obj <- Surv(time = time_vec, event = event_vec)
  fit <- survfit(surv_obj ~ 1)
  
  times <- fit$time
  surv_probs <- fit$surv
  
  cumhaz <- -log(surv_probs)
  
  delta_time <- diff(c(0, times))
  delta_hazard <- diff(c(0, cumhaz))
  hazard_rate <- delta_hazard / delta_time

  df_haz = data.frame(time=unique(time_vec),hazard=hazard_rate,
                      time1= beta_vec[1] * unique(time_vec)^(beta_vec[1] - 1),
                      time3= beta_vec[3] * unique(time_vec)^(beta_vec[3] - 1)
                      )
  df_haz[(nrow(df_haz)),"hazard"] = df_haz[(nrow(df_haz)-1),"hazard"]
  n = nrow(df_haz)  
  subset_df1 = df_haz %>% filter(time<max(times)*ratio1)
  subset_df3 = df_haz %>% filter(time>max(times)*ratio3)
  subset_df2 = df_haz %>% filter(time>max(times)*ratio1 & time<max(times)*ratio3 )
  fit1 <- lm(hazard ~ 0 + time1, data = subset_df1)  # '0 +'는 intercept 제외
  fit3 <- lm(hazard ~ 0 + time3, data = subset_df3)  # '0 +'는 intercept 제외
  
  # 결과 확인
  lambda_est1 <- coef(fit1)[1]
  lambda_est3 <- coef(fit3)[1]
  lambda_est2 = mean(subset_df2[,"hazard"])
  lambda_vec = c(lambda_est1,lambda_est2,lambda_est3)
  return(lambda_vec)
}

