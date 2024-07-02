library(tidyverse)
library(mgcv)
library(splines)
library(survey)


#----- Simulate natural response curve
#===============================================================================
# Given the response probability starts at zero, rises to a peak, and then gradually declines.
# t is set from 0 to 24 hours, the response ranges from 0 to 1
# ==> Weibull distribution or log-normal distribution may be straightforward 
# Here I use Weibull distribution 
# I assume the response peak is at hours 3-5 with about maximum of 75% of participants response

# The scale parameter is set at 15 ($\lambda$ = 15)
# Shape parameter is set at 1.3 ($\kappa$ = 1.3)
# $Y_0 \sim Weibull(\lambda = 15, \kappa = 1.3)$
# A scale factor to make sure the prob reachs a peak of 0.75

# lambda <- 15
# k <- 1.3
# 
# 
# t <- 0:24
# 
# # Weibull PDF
# weibull_pdf <- dweibull(t, shape = k, scale = lambda)
# 
# # Scale the y-axis to have a peak around 0.75
# scaling_factor <- 0.75/max(weibull_pdf)
# weibull_pdf_scaled <- weibull_pdf * scaling_factor
# 
# df_WB <- data.frame(t = t, response = weibull_pdf_scaled)
# 
# # Plot the Weibull distribution
# df_WB |>
#   ggplot(aes(x = t, y = response)) +
#   geom_line(linewidth = 1, color = "red") +
#   theme_minimal() +
#   scale_y_continuous(limit = c(0, 1))




#----- Function 1: function for natural response curve
#===============================================================================
gen_WB_curve <- function(t = 24, peak = 0.75, lambda = 15, k = 1.3) {
  x <- seq(0, t)
  # Weibull PDF
  weibull_pdf <- dweibull(x, shape = k, scale = lambda)
  
  # Scale the y-axis to have a peak around 0.75
  scaling_factor <- peak / max(weibull_pdf)
  weibull_pdf_scaled <- weibull_pdf * scaling_factor
  
  df_WB <- data.frame(t = x, response = weibull_pdf_scaled)
  return(df_WB)
}


#----- Function 2: Simulate data for 3 scenarios
#===============================================================================

genData <- function(n, scen, t = 24, peak = 0.75, lambda = 15, k = 1.3) {
  
  # Generate natural response curve
  df_WB <- gen_WB_curve(t = t, peak = peak, lambda = lambda, k = k)
  u    <- rnorm(n, 0, 2)
  x1   <- rbinom(n, size = 1, prob = plogis(0.5*u))
  x2   <- 0.5*u + rnorm(n, 0, 0.1)
  x3   <- rnorm(n, 0, 2)
  
  # Set up treatment A and outcome Y under different scenarios
  df_mat <- matrix(NA, n, (t+1)*3 + 4)
  
  df_mat[, (t+1)*3 + 2] <- x1
  df_mat[, (t+1)*3 + 3] <- x2
  df_mat[, (t+1)*3 + 4] <- x3
  
  if (scen == "small") {
    A    <- rbinom(n, size=1, prob = plogis(0.5*x1 + 0.5*x3))
    
    df_mat[, (t+1)*3 + 1] <- A
    
    for (i in 2:(t+1)) {
      # counterfactual outcome for A = 1
      yt1 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 1 + 0.5*u))
      
      # counterfactual outcome for A = 0
      yt0 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 0 + 0.5*u))
      
      # Observed outcome
      yt <- yt1*A + yt0*(1 - A)
      
      df_mat[, 3*(i-2) + 4] <- yt1
      df_mat[, 3*(i-2) + 5] <- yt0
      df_mat[, 3*(i-2) + 6] <- yt
      
    }
    df_mat[, 1] <- rep(0, n)
    df_mat[, 2] <- rep(0, n)
    df_mat[, 3] <- rep(0, n)
  } else if (scen == "medium") {
    A    <- rbinom(n, size=1, prob = plogis(0.5*x1 + 0.5*x3))
    
    df_mat[, (t+1)*3 + 1] <- A
    
    for (i in 2:(t+1)) {
      # counterfactual outcome for A = 1
      yt1 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 1 + 0.5*u + 0.5*x3))
      
      # counterfactual outcome for A = 0
      yt0 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 0 + 0.5*u + 0.5*x3))
      
      # Observed outcome
      yt <- yt1*A + yt0*(1 - A)
      
      df_mat[, 3*(i-2) + 4] <- yt1
      df_mat[, 3*(i-2) + 5] <- yt0
      df_mat[, 3*(i-2) + 6] <- yt
      
    }
    df_mat[, 1] <- rep(0, n)
    df_mat[, 2] <- rep(0, n)
    df_mat[, 3] <- rep(0, n)

  } else {
    A    <- rbinom(n, size=1, prob = plogis(0.5*x1 + 0.5*x3 + 0.5*x2))
    
    df_mat[, (t+1)*3 + 1] <- A
    
    for (i in 2:(t+1)) {
      # counterfactual outcome for A = 1
      yt1 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 1 + 0.5*u + 0.75*x3 + 0.75*x2))
      
      # counterfactual outcome for A = 0
      yt0 <- rbinom(n, size=1, prob = plogis(qlogis(df_WB$response[i]) + 0 + 0.5*u + 0.75*x3 + 0.75*x2))
      
      # Observed outcome
      yt <- yt1*A + yt0*(1 - A)
      
      df_mat[, 3*(i-2) + 4] <- yt1
      df_mat[, 3*(i-2) + 5] <- yt0
      df_mat[, 3*(i-2) + 6] <- yt
      
    }
    df_mat[, 1] <- rep(0, n)
    df_mat[, 2] <- rep(0, n)
    df_mat[, 3] <- rep(0, n)
  }
  
  # Return data.frame
  colnames(df_mat) <- c("Y1_0", "Y0_0", "Y_0",
                        paste0(c("Y1_", "Y0_", "Y_"), rep(1:t, each = 3)),
                        "A", "x1", "x2", "x3")
  df_mat <- as.data.frame(df_mat)
  df_mat$id <- 1:n

  # Data for counterfactual outcome of A = 1
  df_Y1 <- df_mat |> select(id, starts_with("Y1_")) |>
    gather(-c(id), key = "t", value = "Y1") |>
    mutate(t = as.numeric(str_remove_all(t, "Y1_"))) |> 
    arrange(id, t)
  
  # Data for counterfactual outcome of A = 0
  df_Y0 <- df_mat |> select(id, starts_with("Y0_")) |>
    gather(-c(id), key = "t", value = "Y0") |>
    mutate(t = as.numeric(str_remove_all(t, "Y0_"))) |> 
    arrange(id, t)
  
  # Data for observed outcome of 
  df_Y <- df_mat |> select(id, A, x1, x2, x3, starts_with("Y_")) |>
    gather(-c(id, A, x1, x2, x3), key = "t", value = "Y") |>
    mutate(t = as.numeric(str_remove_all(t, "Y_"))) |>
    arrange(id, t)
  
  df_all <- df_Y |> 
    left_join(df_Y1, by = c("id", "t")) |>
    left_join(df_Y0, by = c("id", "t"))
  
 return(df_all)
}




#----- Function 3: plot simulated data
#===============================================================================
plot_genData <- function(n = 1000, scen, t = 24, peak = 0.7, lambda = 15, k = 1.3) {
  indat <- genData(n = 1000, scen = scen, t = 24, peak = 0.7, lambda = 15, k = 1.3)
  
  f1_ob <- indat |>
    group_by(t, A) |>
    summarise(prob = mean(Y)) |>
    ggplot(aes(x = as.factor(t), y = prob, color = as.factor(A), group = as.factor(A))) + 
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Time", y = "Response (%)", color = "Group",
         title = "Observed reponse curves (biased)") +
    theme_bw()
  
  A1_sum <- indat |>
    group_by(t) |>
    summarise(prob = mean(Y1), group = "1")
  
  A0_sum <- indat |>
    group_by(t) |>
    summarise(prob = mean(Y0), group = "0")
  
  f1_cft <- rbind(A1_sum, A0_sum) |>
    ggplot(aes(x = as.factor(t), y = prob, color = group, group = group)) + 
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_brewer(palette = "Set1") +
    labs(x = "Time", y = "Response (%)", color = "Group",
         title = "Counterfactual reponse curves (true)") +
    theme_bw()
  return(cowplot::plot_grid(f1_cft, f1_ob, labels = "AUTO", ncol = 2))
}


plot_genData(n = 10000, scen = "small")
plot_genData(n = 10000, scen = "medium")
plot_genData(n = 10000, scen = "large")






#----- Test bias for 3 scenarios (i.e., bias with no adjustment)
#===============================================================================
#----- Simplify the function to calculate area under the curve only (no SE)

auc_smooth_trap_area <- function(a, b, dat, resvar, tvar = "t") {
  indat <- dat
  
  indat$pred <- indat[, resvar]
  indat$t <- indat[, tvar]
  
  indat_sum <- indat |> 
    group_by(t) |>
    summarise(prob = mean(pred)) |>
    as.data.frame()
  
  intev        <- b - a
  idx          <- a <= indat_sum$t & indat_sum$t <= b
  time         <- indat_sum$t[idx]
  pred         <- indat_sum$prob[idx]
  # AUC using trapezoidal rule (pracma::trapz(time[ord],pred[ord])/intev)
  ord          <- order(time)
  area         <- sum(diff(time[ord])*zoo::rollmean(pred[ord], 2))/intev
  return(area)
}



#----- Small confounding
set.seed(12345)
df_small        <- genData(n = 1000, scen = "small", t = 24, peak = 0.7, lambda = 15, k = 1.3)
trueAUC_Y1_s    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_small, resvar = "Y1")
trueAUC_Y0_s    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_small, resvar = "Y0")
(trueAUC_diff_s <- trueAUC_Y1_s - trueAUC_Y0_s)
biasAUC_Y1_s    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_small |> filter(A == 1), resvar = "Y")
biasAUC_Y0_s    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_small |> filter(A == 0), resvar = "Y")
(biasAUC_diff_s <- biasAUC_Y1_s - biasAUC_Y0_s)
biasAUC_diff_s  - trueAUC_diff_s


#----- Medium confounding
df_medium       <- genData(n = 1000, scen = "medium", t = 24, peak = 0.7, lambda = 15, k = 1.3)
trueAUC_Y1_m    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_medium, resvar = "Y1")
trueAUC_Y0_m    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_medium, resvar = "Y0")
(trueAUC_diff_m <- trueAUC_Y1_m - trueAUC_Y0_m)
biasAUC_Y1_m    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_medium |> filter(A == 1), resvar = "Y")
biasAUC_Y0_m    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_medium |> filter(A == 0), resvar = "Y")
(biasAUC_diff_m <- biasAUC_Y1_m - biasAUC_Y0_m)
biasAUC_diff_m - trueAUC_diff_m


#----- Large confounding
df_large <- genData(n = 1000, scen = "large", t = 24, peak = 0.7, lambda = 15, k = 1.3)
trueAUC_Y1_l    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_large, resvar = "Y1")
trueAUC_Y0_l    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_large, resvar = "Y0")
(trueAUC_diff_l <- trueAUC_Y1_l - trueAUC_Y0_l)
biasAUC_Y1_l    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_large |> filter(A == 1), resvar = "Y")
biasAUC_Y0_l    <- auc_smooth_trap_area(a = 0, b = 24, dat = df_large |> filter(A == 0), resvar = "Y")
(biasAUC_diff_l <- biasAUC_Y1_l - biasAUC_Y0_l)
biasAUC_diff_l - trueAUC_diff_l




#----- Function 4: Compare different adjustment methods
#===============================================================================
doSimulation <- function(n = 1000, scen, t = 24, peak = 0.7, lambda = 15, k = 1.3, 
                         a = 0, b = 24, boot = NULL, conf_lv = 0.95, seed = NULL) { 
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  if (scen == "small") {
    dat         <- genData(n = n, scen = "small", t = t, peak = peak, lambda = lambda, k = k)
    formula1    <- Y ~  x1 + s(t)
    formula2    <- A ~  x1 + x3
    formula3    <- Y ~  x1 + bs(t, df = 5) # For svyglm() as mgcv cannot use ipw
  } else if (scen == "medium") {
    dat         <- genData(n = n, scen = "medium", t = t, peak = peak, lambda = lambda, k = k)
    formula1    <- Y ~  x1 + x3 + s(t)
    formula2    <- A ~  x1 + x3
    formula3    <- Y ~  x1 + x3 + bs(t, df = 5) # For svyglm() as mgcv cannot use ipw
  } else {
    dat         <- genData(n = n, scen = "large", t = t, peak = peak, lambda = lambda, k = k)
    formula1    <- Y ~  x1 + x2 + x3 + s(t)
    formula2    <- A ~  x1 + x2 + x3
    formula3    <- Y ~  x1 + x2 + x3 + bs(t, df = 5) # For svyglm() as mgcv cannot use ipw
  }
  
  trueAUC_Y1    <- auc_smooth_trap_area(a = a, b = b, dat = dat, resvar = "Y1")
  trueAUC_Y0    <- auc_smooth_trap_area(a = a, b = b, dat = dat, resvar = "Y0")
  trueAUC_diff  <- trueAUC_Y1 - trueAUC_Y0
  
  biasAUC_Y1    <- auc_smooth_trap_area(a = a, b = b, dat = dat[dat$A==1,], resvar = "Y")
  biasAUC_Y0    <- auc_smooth_trap_area(a = a, b = b, dat = dat[dat$A==0,], resvar = "Y")
  biasAUC_diff  <- biasAUC_Y1 - biasAUC_Y0
  
  
  #----- Parametric G-formula
  par_gform <- function(dat, ...) {
    indat       <- dat
    gam1        <- gam(formula1, data = indat[indat$A==1,], method = "REML", family = binomial) 
    gam0        <- gam(formula1, data = indat[indat$A==0,], method = "REML", family = binomial)
    indat$pred1 <- predict(gam1, type = "response", newdata = indat |> mutate(A = 1))
    indat$pred0 <- predict(gam0, type = "response", newdata = indat |> mutate(A = 0))
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred1")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred0")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  par_gform_auc <- par_gform(dat = dat)
  
  
  #----- IPW
  iptw <- function(dat, ...) {
    indat       <- dat
    ps_mod      <- glm(formula2, data=indat, family="binomial")
    pscore      <- ifelse(indat$A==0, 1-predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
    indat$w     <- 1/pscore
    glm_w       <- svyglm(Y ~ A + bs(t, df = 5), family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w, data = indat))
    indat$pred  <- predict(glm_w, newdata = indat, type = "response")
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat[indat$A==1,], resvar = "pred")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat[indat$A==0,], resvar = "pred")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  
  iptw_auc <- iptw(dat = dat)
  
  #----- Double-robust methods
  DB_both <- function(dat,...) {
    indat       <- dat
    ps_mod      <- glm(formula2, data=indat, family="binomial")
    pscore      <- ifelse(indat$A==0, 1-predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
    indat$w     <- 1/pscore
    glm1        <- svyglm(formula3, family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==1], data = indat[indat$A==1,]))
    glm0        <- svyglm(formula3, family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==0], data = indat[indat$A==0,]))
    indat$pred1 <- predict(glm1, type = "response", newdata = indat |> mutate(A = 1))
    indat$pred0 <- predict(glm0, type = "response", newdata = indat |> mutate(A = 0))
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred1")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred0")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  
  DB_auc <- DB_both(dat = dat)
  
  
  #----- Double-robust methods: Misclassification of treatment model
  DB_A <- function(dat,...) {
    indat       <- dat
    ps_mod      <- glm(A ~ 1, data=indat, family="binomial")
    pscore      <- ifelse(indat$A==0, 1-predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
    indat$w     <- 1/pscore
    glm1        <- svyglm(formula3, family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==1], data = indat[indat$A==1,]))
    glm0        <- svyglm(formula3, family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==0], data = indat[indat$A==0,]))
    indat$pred1 <- predict(glm1, type = "response", newdata = indat |> mutate(A = 1))
    indat$pred0 <- predict(glm0, type = "response", newdata = indat |> mutate(A = 0))
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred1")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred0")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  
  DB_A_auc <- DB_A(dat = dat)
  
  #----- Double-robust methods: Misclassification of outcome model
  DB_Y <- function(dat,...) {
    indat       <- dat
    ps_mod      <- glm(formula2, data=indat, family="binomial")
    pscore      <- ifelse(indat$A==0, 1-predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
    indat$w     <- 1/pscore
    glm1        <- svyglm(Y ~ bs(t, df = 5), family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==1], data = indat[indat$A==1,]))
    glm0        <- svyglm(Y ~ bs(t, df = 5), family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==0], data = indat[indat$A==0,]))
    indat$pred1 <- predict(glm1, type = "response", newdata = indat |> mutate(A = 1))
    indat$pred0 <- predict(glm0, type = "response", newdata = indat |> mutate(A = 0))
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred1")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred0")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  
  DB_Y_auc <- DB_Y(dat = dat)
  
  #----- Double-robust methods: Misclassification of both treatment and outcome models
  DB_AY <- function(dat,...) {
    indat       <- dat
    ps_mod      <- glm(A ~ 1, data=indat, family="binomial")
    pscore      <- ifelse(indat$A==0, 1-predict(ps_mod, type = "response"), predict(ps_mod, type = "response"))
    indat$w     <- 1/pscore
    glm1        <- svyglm(Y ~ bs(t, df = 5), family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==1], data = indat[indat$A==1,]))
    glm0        <- svyglm(Y ~ bs(t, df = 5), family = "quasibinomial", 
                          design = svydesign(~ 1, weights = ~ indat$w[indat$A==0], data = indat[indat$A==0,]))
    indat$pred1 <- predict(glm1, type = "response", newdata = indat |> mutate(A = 1))
    indat$pred0 <- predict(glm0, type = "response", newdata = indat |> mutate(A = 0))
    adj_AUC1    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred1")
    adj_AUC0    <- auc_smooth_trap_area(a = a, b = b, dat = indat, resvar = "pred0")
    AUC_diff    <- adj_AUC1 - adj_AUC0
    return(list(adj_AUC1 = adj_AUC1, adj_AUC0 = adj_AUC0, AUC_diff = AUC_diff))
  }
  
  DB_AY_auc <- DB_AY(dat = dat)
  
  
  df_est <- data.frame(
    method = c("Gform", "IPW", "DB", "DB_A", "DB_Y", "DB_AY"),
    trueAUC_Y1 = trueAUC_Y1,
    trueAUC_Y0 = trueAUC_Y0,
    biasAUC_Y1 = biasAUC_Y1,
    biasAUC_Y0 = biasAUC_Y0,
    adj_AUC1 = c(par_gform_auc$adj_AUC1, iptw_auc$adj_AUC1, DB_auc$adj_AUC1, DB_A_auc$adj_AUC1, DB_Y_auc$adj_AUC1, DB_AY_auc$adj_AUC1),
    adj_AUC0 = c(par_gform_auc$adj_AUC0, iptw_auc$adj_AUC0, DB_auc$adj_AUC0, DB_A_auc$adj_AUC0, DB_Y_auc$adj_AUC0, DB_AY_auc$adj_AUC0),
    true_AUC_diff = trueAUC_diff,
    bias_AUC_diff = biasAUC_diff,
    adj_AUC_diff = c(par_gform_auc$AUC_diff, iptw_auc$AUC_diff, DB_auc$AUC_diff, DB_A_auc$AUC_diff, DB_Y_auc$AUC_diff, DB_AY_auc$AUC_diff))
  
  if (!is.null(boot)) {
    
    #----- Bootstrap
    #=========================================================================
    df_result_b <- NULL
    
    
    #--- Bootstrap data
    for (i in 1:boot) {
      idx      <- sample(1:n, n, replace = T)
      freq_id  <- table(idx)
      df_boot  <- NULL
      
      for(j in 1:max(freq_id)) {
        # Loop over repeated id
        temp_df <- dat[dat$id %in% names(freq_id[freq_id %in% c(j:max(freq_id))]), ]
        temp_df$boot_id <- paste0(temp_df$id, "_", j)
        df_boot <- rbind(df_boot, temp_df)
      }
      
      trueAUC_Y1_b    <- auc_smooth_trap_area(a = a, b = b, dat = df_boot, resvar = "Y1")
      trueAUC_Y0_b    <- auc_smooth_trap_area(a = a, b = b, dat = df_boot, resvar = "Y0")
      trueAUC_diff_b  <- trueAUC_Y1_b - trueAUC_Y0_b
      
      biasAUC_Y1_b    <- auc_smooth_trap_area(a = a, b = b, dat = df_boot[df_boot$A==1,], resvar = "Y")
      biasAUC_Y0_b    <- auc_smooth_trap_area(a = a, b = b, dat = df_boot[df_boot$A==0,], resvar = "Y")
      biasAUC_diff_b  <- biasAUC_Y1_b - biasAUC_Y0_b
      
      
      #--- Bootstrap estimates
      
      # parametric G-formula
      par_gform_b <- par_gform(dat = df_boot)
      
      # IPW
      iptw_b <- iptw(dat = df_boot)
      
      # DB both
      DB_b <- DB_both(dat = df_boot)
      
      # Double-robust methods: Misclassification of treatment model 
      DB_A_b <- DB_A(dat = df_boot)
      
      # Double-robust methods: Misclassification of outcome model
      DB_Y_b <- DB_Y(dat = df_boot)
      
      # Double-robust methods: Misclassification of both treatment and outcome models
      DB_AY_b <- DB_AY(dat = df_boot)
      
      df_b <- data.frame(
        method = c("Gform", "IPW", "DB", "DB_A", "DB_Y", "DB_AY"),
        trueAUC_Y1 = trueAUC_Y1_b,
        trueAUC_Y0 = trueAUC_Y0_b,
        biasAUC_Y1 = biasAUC_Y1_b,
        biasAUC_Y0 = biasAUC_Y0_b,
        adj_AUC1 = c(par_gform_b$adj_AUC1, iptw_b$adj_AUC1, DB_b$adj_AUC1, DB_A_b$adj_AUC1, DB_Y_b$adj_AUC1, DB_AY_b$adj_AUC1),
        adj_AUC0 = c(par_gform_b$adj_AUC0, iptw_b$adj_AUC0, DB_b$adj_AUC0, DB_A_b$adj_AUC0, DB_Y_b$adj_AUC0, DB_AY_b$adj_AUC0),
        true_AUC_diff = trueAUC_diff_b,
        bias_AUC_diff = biasAUC_diff_b,
        adj_AUC_diff = c(par_gform_b$AUC_diff, iptw_b$AUC_diff, DB_b$AUC_diff, DB_A_b$AUC_diff, DB_Y_b$AUC_diff, DB_AY_b$AUC_diff),
        iteration = i)
      
      df_result_b <- rbind(df_result_b, df_b)
    }
    
    # Summary data for bootstrap
    df_est <- df_result_b |> 
      group_by(method) |>
      summarise(trueAUC_Y1_lb     = quantile(trueAUC_Y1, probs = (1 - conf_lv)/2),
                trueAUC_Y1_ub     = quantile(trueAUC_Y1, probs = (1 + conf_lv)/2),
                trueAUC_Y0_lb     = quantile(trueAUC_Y0, probs = (1 - conf_lv)/2),
                trueAUC_Y0_ub     = quantile(trueAUC_Y0, probs = (1 + conf_lv)/2),
                biasAUC_Y1_lb     = quantile(biasAUC_Y1, probs = (1 - conf_lv)/2),
                biasAUC_Y1_ub     = quantile(biasAUC_Y1, probs = (1 + conf_lv)/2),
                biasAUC_Y0_lb     = quantile(biasAUC_Y0, probs = (1 - conf_lv)/2),
                biasAUC_Y0_ub     = quantile(biasAUC_Y0, probs = (1 + conf_lv)/2),
                adj_AUC1_lb       = quantile(adj_AUC1, probs = (1 - conf_lv)/2),
                adj_AUC1_ub       = quantile(adj_AUC1, probs = (1 + conf_lv)/2),
                adj_AUC0_lb       = quantile(adj_AUC0, probs = (1 - conf_lv)/2),
                adj_AUC0_ub       = quantile(adj_AUC0, probs = (1 + conf_lv)/2),
                true_AUC_diff_lb  = quantile(true_AUC_diff, probs = (1 - conf_lv)/2),
                true_AUC_diff_ub  = quantile(true_AUC_diff, probs = (1 + conf_lv)/2),
                bias_AUC_diff_lb  = quantile(bias_AUC_diff, probs = (1 - conf_lv)/2),
                bias_AUC_diff_ub  = quantile(bias_AUC_diff, probs = (1 + conf_lv)/2),
                adj_AUC_diff_lb   = quantile(adj_AUC_diff, probs = (1 - conf_lv)/2),
                adj_AUC_diff_ub   = quantile(adj_AUC_diff, probs = (1 + conf_lv)/2))|> 
      ungroup() |>
      right_join(df_est, by = "method") |>
      datawizard::data_rotate(colnames = TRUE, rownames = "param")
  }
  
  return(df_est)
  
}


#----- Estimates
doSimulation(scen = "small", seed = 12345)
doSimulation(scen = "medium", seed = 12345)
doSimulation(scen = "large", seed = 12345)


#----- Estimates with Bootstrap CIs
doSimulation(scen = "small", boot = 5000, seed = 12345)
doSimulation(scen = "medium", boot = 5000, seed = 12345)
doSimulation(scen = "large", boot = 5000, seed = 12345)




#----- Run the loop to calculate average bias
#===============================================================================
ave_bias_tab <- NULL

for (k in c("small", "medium", "large")) {
  for (i in 1: 10) {
    df_out          <- doSimulation(n = 1000, scen = k, seed = i)
    df_out$scenario <- k
    df_out$iter     <- i
    ave_bias_tab    <- rbind(ave_bias_tab, df_out)
  }
}


ave_bias_tab |> 
  group_by(scenario, method) %>%
  summarise(trueAUC_Y1 = mean(trueAUC_Y1),
            trueAUC_Y0 = mean(trueAUC_Y0),
            biasAUC_Y1 = mean(biasAUC_Y1),
            biasAUC_Y0  = mean(biasAUC_Y0),
            adj_AUC1  = mean(adj_AUC1),
            adj_AUC0 = mean(adj_AUC0),
            true_AUC_diff = mean(true_AUC_diff),
            bias_AUC_diff = mean(bias_AUC_diff),
            adj_AUC_diff = mean(adj_AUC_diff))






