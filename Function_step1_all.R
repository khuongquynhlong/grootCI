


#----- Function 0: set up natural response curve
#===============================================================================
dat_guide_funct <- function(n_measure = 25, arm_name = c("Treatment", "Control"), 
                            diff_ratio = 0.08) {
  
  # Data to guide simulation
  xdat <- c(0, 1, 3, 5, 7, 9, 12, 15, 18, 21, 24)
  ydat <- c(9, 3, 2, 1, 1, 1.2, 1.5, 1.8, 2, 2.5, 3)
  
  data <- data.frame(x = xdat, y = ydat)
  
  # Fit a loess
  fit <- loess(y ~ x, data = data)
  
  # Generate data for plotting the fitted curve
  x_fit <- 1:(n_measure-1)
  
  # Mean group A
  y_mA <- predict(fit, newdata = data.frame(x = x_fit))
  y_mB <- y_mA + y_mA*diff_ratio*(abs(x_fit - mean(x_fit)+5)) + rnorm(length(y_mA), 0, 0.1)
  
  # Add baseline value
  y_mA <- c(9, y_mA)
  y_mB <- c(9, y_mB)
  
  df <- data.frame(arm = rep(arm_name, each = n_measure), 
                   y = c(y_mA, y_mB), 
                   x = rep(0:(n_measure-1), 2))
  
  my_plot <- df |> 
    ggplot(aes(x = as.factor(x), y = y, color = arm, group = arm)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(breaks = 1:10) +
    labs(x = "time", y = "NRS score") +
    theme_bw()
  
  return(list(y_mA = y_mA, y_mB = y_mB, plot = my_plot))
}

# x <- dat_guide_funct(n_measure = 25, diff_ratio = 0.1, arm_name = c("Treatment", "Control"))
# x$plot





#----- Function 1a: to Simulate data (wide data)
#===============================================================================
sim_dat_wide_nocov <- function(N = 100, n_measure = 25, arm_name = c("A", "B"), diff_ratio = 0.08) {
  dat_guide <- dat_guide_funct(n_measure = n_measure, arm_name = arm_name, diff_ratio = diff_ratio)
  
  y_mA <- dat_guide$y_mA
  y_mB <- dat_guide$y_mB
  
  df_matA <- df_matB <- matrix(ncol = n_measure, nrow = N)
  
  for (i in 1:n_measure) {
    df_matA[, i] <- extraDistr::rtpois(n = N, lambda = y_mA[i], b = 10)
    df_matB[, i] <- extraDistr::rtpois(n = N, lambda = y_mB[i], b = 10)
  }
  
  df_mat <- rbind(df_matA, df_matB) |> as.data.frame()
  names(df_mat) <- paste0("t", 0:(n_measure - 1))
  df_mat$id <- 1:(2*N)
  df_mat$group <- rep(arm_name, each = N)
  
  return(df_mat)
}

# df_wide_sim <- sim_dat_wide_nocov(N = 1000, n_measure = 25, diff_ratio = 0.08, arm_name = c("Treatment", "Control"))
# rbind(head(df_wide_sim), tail(df_wide_sim))


#----- Function 1b: Create response dataset (change from baseline i.e., ≥50% improvement)
#===============================================================================
sim_dat_change_long_nocov <- function(cutoff = -0.5, dat) 
{
  # generate data
  df <- dat
  
  for (i in 1:nrow(df)) {
    for (j in 2:(ncol(df)-2)) {
      if (df[i, 1] > 0) {
        df[i, j] <- ifelse((df[i, j] - df[i, 1])/df[i, 1] <= cutoff, 1, 0)
      } else {
        df[i, j] <- ifelse((df[i, j] == df[i, 1]), 1, 0)
      }
    }
  }
  
  df[, 1] <- 0
  
  df_change_long <- df |>
    gather(-c(id, group), key = "t", value = "prob") |>
    mutate(t = parse_number(t)) |>
    ungroup() |> as.data.frame()
  
  return(df_change_long)
}

# df_change <- sim_dat_change_long_nocov(cutoff = -0.5, dat = df_wide_sim)
# rbind(head(df_change_long_sim), tail(df_change_long_sim))


#----- Function 3: THE GROUP RESPONSE OUTCOME (GRO) and Area under the curve
#===============================================================================

plot_GRO <- function(dat) {
  indat <- dat
  plot_step <- indat |> 
    group_by(group, t) |>
    summarise(prob = mean(prob)) |> ungroup() |>
    ggplot(aes(x = as.factor(t), y = prob, color = group, group = group)) +
    geom_step(linewidth = 1, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    guides(col = guide_legend(override.aes = list(alpha = 1))) + 
    labs(x = "time", y = "Response (%)") +
    theme_bw()
  plot_line <- indat |> 
    group_by(group, t) |>
    summarise(prob = mean(prob)) |> ungroup() |>
    ggplot(aes(x = as.factor(t), y = prob, color = group, group = group)) +
    geom_line(linewidth = 1, alpha = 0.7) +
    scale_color_brewer(palette = "Set1") +
    guides(col = guide_legend(override.aes = list(alpha = 1))) + 
    labs(x = "time", y = "Response (%)") +
    theme_bw()
  return(list(step = plot_step, line = plot_line))
}

# x <- plot_GRO(dat = df_change)
# x$step
# x$line





#----- Function 3a: GROUP RESPONSE OUTCOME OVER TIME (GROOT) based on GRO (i.e., average curve)
#===============================================================================

auc_funct <- function(a = 1, b, dat, idvar = "id", respvar = "prob", tvar = "t", grvar = "group") {
  
  auc_nonpar <- function(a = 1, b, dat, arm) {
    dat$group <- dat[, grvar]
    
    indat <- dat[dat$group == arm, ]
    
    indat$id <- indat[, idvar]
    indat$prob <- indat[, respvar]
    indat$t <- indat[, tvar]
    
    indat_sum <- indat |> 
      group_by(t) |>
      summarise(prob = mean(prob, na.rm = T),
                n = n(),
                var = (prob*(1-prob))/(n-1)) |>
      as.data.frame()
    
    # Calculate AUC as step function
    intev       <- b - a
    idx         <- a <= indat_sum$t & indat_sum$t <= b
    time        <- c(indat_sum$t[idx], b)
    resp        <- indat_sum$prob[idx]
    var_id      <- indat_sum$var[idx]
    time_diff   <- diff(c(time))
    areas_diff  <- time_diff * resp/intev
    area        <- sum(areas_diff)
    # Check!!
    var_area    <- sum((time_diff^2)*var_id)/(intev^2)
    return(list(area = area, var = var_area))
  }
  
  
  multiarm <- function(val){
    auc_nonpar(arm   = val,
               a     = a,
               b     = b,
               dat   = dat)
  }
  
  
  arm_name <- as.character(unique(sort(dat$group)))
  area_mat <- t(sapply(arm_name, multiarm))
  
  # Set up result matrix
  res_mat <- cbind(area_mat, matrix(NA, length(arm_name), 6))
  res_mat[, 3] <- unlist(res_mat[, 1]) - sqrt(unlist(res_mat[, 2]))*qnorm(1 - 0.05/2)
  res_mat[, 4] <- unlist(res_mat[, 1]) + sqrt(unlist(res_mat[, 2]))*qnorm(1 - 0.05/2)
  
  if (length(arm_name) > 1) {
    for (j in 2:length(arm_name)) {
      # Difference
      res_mat[j, 5] <- area_mat[[1, 1]] - area_mat[[j, 1]]
      # SE difference
      res_mat[j, 6] <- sqrt(area_mat[[1, 2]] + area_mat[[j, 2]])
      # z-score
      res_mat[j, 7] <- res_mat[[j, 5]]/res_mat[[j, 6]]
      # p-value (1 - pchisq(res_mat[[2, 7]]^2, 1))
      res_mat[j, 8] <- pnorm(-abs(res_mat[[j, 7]])) * 2
    }
  }
  
  colnames(res_mat) <- c("area", "var", "lb", "ub", "diff", "se_diff", "z", "pval")
  
  return(res_mat)
}

# auc_res <- auc_funct(a = 1, b = 25, dat = df_change)
# auc_res






#----- Function 3b: GROOT based on Participant-Level Outcome (i.e., Mean ROOT)
#===============================================================================

auc_funct_idv <- function(a = 1, b, dat, idvar = "id", respvar = "prob", tvar = "t", grvar = "group") {
  
  auc_nonpar <- function(a = 1, b, dat, arm) {
    dat$group <- dat[, grvar]
    
    indat <- dat[dat$group == arm, ]
    
    indat$id <- indat[, idvar]
    indat$prob <- indat[, respvar]
    indat$t <- indat[, tvar]
    
    # Define lag variable to indicate time-interval
    intev       <- b - a
    
    
    indat_sum <- indat |> 
      arrange(id, t) |>
      group_by(id) |>
      mutate(tlag = lag(t, default = 0),
             time_diff = t - tlag,
             areas_diff = prob*time_diff) |>
      filter(a <= t & t < b) |> 
      summarise(area_ind = sum(areas_diff)/intev,
                n = n(),
                var_ind = area_ind*(1-area_ind)/(n-1))
    
    
    area        <- mean(indat_sum$area_ind)
    # Check!!!!
    var_area    <- sum(indat_sum$var_ind)/nrow(indat_sum)^2
    
    return(list(area = area, var = var_area))
  }
  
  multiarm <- function(val){
    auc_nonpar(arm   = val,
               a     = a,
               b     = b,
               dat   = dat)
  }
  
  arm_name <- as.character(unique(sort(dat$group)))
  area_mat <- t(sapply(arm_name, multiarm))
  
  # Set up result matrix
  res_mat <- cbind(area_mat, matrix(NA, 2, 6))
  res_mat[, 3] <- unlist(res_mat[, 1]) - sqrt(unlist(res_mat[, 2]))*qnorm(1 - 0.05/2)
  res_mat[, 4] <- unlist(res_mat[, 1]) + sqrt(unlist(res_mat[, 2]))*qnorm(1 - 0.05/2)
  
  if (length(arm_name) > 1) {
    for (j in 2:length(arm_name)) {
      # Difference
      res_mat[j, 5] <- area_mat[[1, 1]] - area_mat[[j, 1]]
      # SE difference
      res_mat[j, 6] <- sqrt(area_mat[[1, 2]] + area_mat[[j, 2]])
      # z-score
      res_mat[j, 7] <- res_mat[[j, 5]]/res_mat[[j, 6]]
      # p-value (1 - pchisq(res_mat[[2, 7]]^2, 1))
      res_mat[j, 8] <- pnorm(-abs(res_mat[[j, 7]])) * 2
    }
  }
  
  colnames(res_mat) <- c("area", "var", "lb", "ub", "diff", "se_diff", "z", "pval")
  
  return(res_mat)
}


# auc_res2 <- auc_funct_idv(a = 1, b = 25, dat = df_change)
# auc_res
# auc_res2



# Check with boostrap


# auc_bootA <- vector(length = 1000)
# auc_bootB <- vector(length = 1000)
# auc_bootDiff <- vector(length = 1000)
# 
# for (i in 1:1000) {
#   idx <- sample(1:2000, 2000, replace = T)
#   bb <- table(idx)
#   
#   boot <- NULL
#   
#   for(zzz in 1:max(bb)) {
#     # Loop over repeated id
#     cc <- df_change[df_change$id %in% names(bb[bb %in% c(zzz:max(bb))]), ]
#     cc$bid <- paste0(cc$id, zzz)
#     boot <- rbind(boot, cc)
#   }
#   
#   auc_res <- auc_funct(a = 1, b = 25, dat = boot)
#   auc_bootA[i] <- auc_res[[1, 1]]
#   auc_bootB[i] <- auc_res[[2, 1]]
#   
#   auc_bootDiff[i] <- auc_bootA[i] - auc_bootB[i]
# }


# mean(auc_bootA)
# var(auc_bootA)
# quantile(auc_bootA, probs = c(0.025, 0.975))
# 
# 
# mean(auc_bootB)
# var(auc_bootB)
# quantile(auc_bootB, probs = c(0.025, 0.975))


