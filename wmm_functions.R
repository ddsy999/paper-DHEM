library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(patchwork)
library(purrr)

############################
# Help functions
############################

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


############################
# Algorithm 
############################

wmm_EM <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   tol = 1e-6,verbose=FALSE
                  ){
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", length(maxGEMiter))
  
  bw = 0
  for (it in 1:maxGEMiter) {
    ### M-step ###
    new_pi = colSums(gamma)/N
    
    new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
    new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = 1
    new_beta = c(new_beta1,new_beta2,new_beta3)
    
    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
  
    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    beta=new_beta; pi=new_pi; lambda=new_lambda;

    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
    ### Stopping rule ###
    if(parameter_diff<tol || it==maxGEMiter){
      if(verbose) cat("EM ","[",it,"]"," beta :",beta ,"\n")
      break
    }
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=1)

    trace[[it]] <- list(pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
  }

  return(list(trace = trace, lambda = lambda, beta = beta))
}


wmm_DAEM <- function(df,
                     theta,
                     maxGEMiter = 1e+3,
                     nsteps = 50,
                     r_init = 0.1,
                     r_end = 1,
                     tol=1e-6,verbose=FALSE
                    ){
  method = "DAEM"
  r_grid <- exp(seq(log(r_init), log(r_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  bw = 0
  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]
    
    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }

  return(list(trace = trace, pi = pi, lambda = lambda, beta = beta))
}

wmm_BM <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   nsteps = 50,
                   bw_init = 1e-1,
                   bw_end = 1e-5,
                   tol=1e-6,verbose=FALSE
                  ){
  method="BM"         
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  r=1

  for( hyperIter in 1:nsteps){

    bw=bw_grid[hyperIter]

    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM = function(df,
                    theta,
                    maxGEMiter=1e+3,
                    nsteps=1e+2,
                    r_init=0.1,
                    r_end=1,
                    bw_init=1e-1,
                    bw_end=1e-5,
                    tol=1e-6,verbose=FALSE
  ){
  method = "DHEM"
  r_grid  <- exp(seq(log(r_init), log(r_end), length.out = nsteps))
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]
    bw=bw_grid[hyperIter]

    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2+(pi-new_pi)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM_adaptive <- function(df,
                              theta,
                              maxGEMiter=1e+3,
                              nsteps=1e+2,
                              r_init=0.1,
                              r_end=1,
                              bw_init=1e-1,
                              eta=0.1,
                              tol=1e-6,verbose=FALSE
                            ){
  method = "adapDHEM"
  # r schedule (r -> 1)
  r_grid <- exp(seq(log(r_init), 0, length.out = nsteps))

  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)

  N=nrow(df)
  K=ncol(gamma)

  time_vec  = df$time
  event_vec = df$event

  trace <- vector("list", nsteps)

  bw = bw_init
  for( hyperIter in 1:nsteps){
    r <- r_grid[hyperIter]
    #bw = bw*0.1
    for( gemIter in 1:maxGEMiter){
    ### M-step ###
    new_pi = colSums(gamma)/N
    
    new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
    new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = 1
    new_beta = c(new_beta1,new_beta2,new_beta3)
    
    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
  
    theta0 = list(pi=pi,beta=beta,lambda=lambda)
    theta1 = list(pi=new_pi,beta=new_beta,lambda=new_lambda)

    # Checking ACC
    acc1 = FALSE;acc2 = FALSE;
    deltaDKL = wmm_delta_DKL(df,theta0,theta1,r)
    DKL      = wmm_DKL(df,theta0,theta1)
    deltaB   = wmm_bar_diff(theta1,theta0)
    #cat(deltaDKL-bw*wmm_bar_diff(theta1,theta0),"\n")
    if(deltaDKL-bw*deltaB<0){
      # Acc 1st test      
      if(!is.finite(deltaDKL)||!is.finite(DKL)||DKL<0) break
      if(deltaDKL<eta*DKL) {
        #cat(deltaDKL,eta*DKL,"\n")
        break}else{acc1=TRUE}
      # Acc 2nd test
      if(bw*abs(deltaB)>eta*DKL){
        bw = min(bw,eta*DKL/abs(deltaB))
        next
      }else{acc2 = TRUE}
    }else{
      acc1=TRUE;acc2=TRUE;
    }
    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    #print(parameter_diff)
    beta=new_beta; pi=new_pi; lambda=new_lambda;
    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
     
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    ### Stopping rule ###
    if(parameter_diff<tol || gemIter==maxGEMiter){
      if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3," para diff: ", parameter_diff,"\n")
      break
    }
    }

    if(acc1&&acc2){
      trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
      last_acc = list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
    }
  }

  return(list(trace = trace,pi = last_acc$pi, lambda = last_acc$lambda, beta = last_acc$beta,bw=last_acc$bw,dQbeta1 = last_acc$dQbeta1, dQbeta3=last_acc$dQbeta3))
}

make_result_df <- function(fit, method, K = 3) {
  non_null_idx  <- which(!sapply(fit$trace, is.null))
  valid_traces  <- fit$trace[non_null_idx]

  data.frame(
    method  = method,
    r       = sapply(valid_traces, function(x) if (!is.null(x$r))  x$r  else 1),
    bw      = sapply(valid_traces, function(x) if (!is.null(x$bw)) x$bw else 0),
    pi1     = sapply(valid_traces, function(x) x$pi[1]),
    pi2     = sapply(valid_traces, function(x) x$pi[2]),
    pi3     = sapply(valid_traces, function(x) x$pi[3]),
    beta1   = sapply(valid_traces, function(x) x$beta[1]),
    beta2   = sapply(valid_traces, function(x) x$beta[2]),
    beta3   = sapply(valid_traces, function(x) x$beta[3]),
    lambda1 = sapply(valid_traces, function(x) x$lambda[1]),
    lambda2 = sapply(valid_traces, function(x) x$lambda[2]),
    lambda3 = sapply(valid_traces, function(x) x$lambda[3]),
    dQbeta1 = sapply(valid_traces, function(x) if (length(x$dQbeta1) == 1) x$dQbeta1 else NA_real_),
    dQbeta3 = sapply(valid_traces, function(x) if (length(x$dQbeta3) == 1) x$dQbeta3 else NA_real_)
  )
}

plot_beta_trace <- function(df, title = NULL,
                            axis_text_y_size  = 14,
                            axis_title_y_size = 14,
                            axis_text_x_size  = 14,
                            axis_title_x_size = 14,
                            title_size        = 20
                          ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,
      angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = r, y = beta1)) +
    geom_line() +
    coord_cartesian(ylim = c(0, 1), xlim = c(r_init, 1)) +
    labs(x = "Annealing parameter", y = expression(beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = r, y = beta3)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1)) +
    labs(x = "Annealing parameter", y = expression(beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


plot_dQbeta_trace <- function(df, title = NULL,
                              axis_text_y_size  = 14,
                              axis_title_y_size = 14,
                              axis_text_x_size  = 14,
                              axis_title_x_size = 14,
                              title_size        = 20
                            ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle =0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = r, y = dQbeta1)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1))+
    labs(x = "Annealing parameter", y = expression(nabla~Q~beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = r, y = dQbeta3)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1))+
    labs(x = "Annealing parameter", y = expression(nabla~Q~beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


best_row <- function(df) {
  idx <- which.min(abs(df$dQbeta1) + abs(df$dQbeta3))
  df[idx, c("method","pi1","pi2","pi3","beta1","beta3","lambda1","lambda2","lambda3","dQbeta1","dQbeta3")]
}


plot_beta_trace_barrier <- function(df, title = NULL,
                            axis_text_y_size  = 14,
                            axis_title_y_size = 14,
                            axis_text_x_size  = 10,
                            axis_title_x_size = 14,
                            title_size        = 20
                          ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,
      angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = bw, y = beta1)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, 1)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, 1)
)+
    labs(x = "Barrier parameter", y = expression(beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = bw, y = beta3)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, 1)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, 1)
)+
    labs(x = "Barrier parameter", y = expression(beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}


plot_dQbeta_trace_barrier <- function(df, title = NULL,
                              axis_text_y_size  = 14,
                              axis_title_y_size = 14,
                              axis_text_x_size  = 10,
                              axis_title_x_size = 14,
                              title_size        = 20
                            ) {

  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle =0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )

  p1 <- ggplot(df, aes(x = bw, y = dQbeta1)) +
    geom_line() +
    # scale_x_reverse() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, 1)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, 1)
)+
    labs(x = "Barrier parameter", y = expression(nabla~Q~beta[1])) +
    base_theme

  p3 <- ggplot(df, aes(x = bw, y = dQbeta3)) +
    geom_line() +
    scale_x_continuous(
  trans = scales::trans_new(
    name = "revlog10",
    transform = function(x) -log10(x),
    inverse   = function(x) 10^(-x)
  ),
  breaks = scales::log_breaks(base = 10)(c(1e-8, 1)),
  labels = scales::label_math(10^.x),
  limits = c(1e-8, 1)
)+
    labs(x = "Barrier parameter", y = expression(nabla~Q~beta[3])) +
    base_theme

  (p1 / p3) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
}

