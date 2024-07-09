library(plyr)
library(tidyverse)
library(haven)
library(mgcv)
library(splines)
library(survey)
library(lme4)
library(lmerTest)
library(doParallel)

# set up for paralell running
cores <- detectCores()
cores
registerDoParallel(cores=cores-4)

#-------------------------------------------------------------------------------
ttrial <- read_dta("Data/Mainpaper.dta")

ttrial <- ttrial |> select(pid, monthvisit, phy_yes, sex_yes, vit_yes, names(ttrial))
ttrial <- ttrial |> mutate(gr = ifelse(arm == "A", "Testosterone", "Placebo"),
                           monthvisit_f = as.factor(monthvisit),
                           tr = ifelse(arm == "A", 1, 0)) 
# For all participants
ttrial <- ttrial |> mutate(all_part = 1)

glimpse(ttrial)

ttrial |> group_by(pid) |> slice(1) |> pull(arm) |> table()


#---------- Check sample (All ok)
#===============================================================================
#--- Men enrolled in Sexual Function Trial

# PDQ-Q4 score
ttrial |> filter(sex_yes == 1 & !is.na(Activity_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

ttrial |> filter(!is.na(chg_Activity)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()


# DISF-M-II sexual desire score
ttrial |> filter(sex_yes == 1 & !is.na(sexdes_der_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# IIEF erectile function score**
ttrial |> filter(sex_yes == 1 & !is.na(IIEFEF_DER_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()


#--- Men enrolled in Physical Function Trial

# increase of ≥50 m in 6-min walk test 
ttrial |> filter(phy_yes == 1 & !is.na(WLK_Distance_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

ttrial |> filter(!is.na(WLK_Distance_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# PF-10 score
ttrial |> filter(phy_yes == 1 & !is.na(pf10_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

ttrial |> filter(!is.na(pf10_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()


#--- Men enrolled in Vitality Trial

# Increase of ≥4 in FACIT– Fatigue score 
ttrial |> filter(vit_yes == 1 & !is.na(FACIT_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

ttrial |> filter(!is.na(FACIT_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# SF-36 vitality score
ttrial |> filter(vit_yes == 1 & !is.na(SF36_IVR_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# PANAS positive affect score
ttrial |> filter(vit_yes == 1 & !is.na(PANASpos_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# PANAS negative affect score
ttrial |> filter(vit_yes == 1 & !is.na(PANASneg_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()

# PHQ-9 depression score
ttrial |> filter(vit_yes == 1 & !is.na(PHQ9_bsl)) |> group_by(pid) |>
  slice(1) |> pull(arm) |> table()


#----- balancing factors
# baseline total testosterone level (≤200 or >200 ng per deciliter): 
table(ttrial$Tunder)
# age (≤75 or >75 years)
table(ttrial$ageunder)
# trial site
table(ttrial$sitenum)
# participation in the main trials
table(ttrial$Tunder)
# use or nonuse of antidepressants
table(ttrial$dep_treat)
# and use or nonuse of phosphodiesterase type 5 inhibitors
table(ttrial$pde_yes)







#---------- Prepare all functions
#===============================================================================

#----------  Set up theme
mytheme <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14,color = "grey10",  face = "bold", hjust = 0.5),
      plot.subtitle = element_text(face = "italic", color = "gray10", size = 14),
      plot.caption = element_text(face = "italic", size = 14, color = "gray10"),
      axis.line = element_line(linetype = "solid"),
      axis.text.x = element_text(color = "gray10", size = 14),
      axis.text.y = element_text(color = "gray10", size = 14),
      # axis.ticks = element_blank(),
      axis.title.x = element_text(color = "gray10", size = 14),
      axis.title.y = element_text(color = "gray10", size = 14),
      # panel.grid.minor = element_blank(),
      # panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 14, face = "bold"),
      legend.direction = "horizontal",
      legend.position = "top",
      legend.background = element_rect(fill = NA, color = NA),
      legend.text = element_text(size = 14),
      legend.key.width = unit(2, "line"),
      strip.text = element_text(size = 14, face = "bold"),
      strip.background = element_rect(fill = NA, color = NA)
    )
}

#===============================================================================
#----- Simplify the function to calculate area under the curve only (no SE)
auc_smooth_trap_ttrial <- function(a, b, dat, resvar, tvar = "t") {
  indat <- as.data.frame(dat)
  
  indat$pred <- indat[, resvar]
  indat$t <- indat[, tvar]
  
  indat_sum <- indat |> 
    group_by(t) |>
    summarise(prob = mean(pred, na.rm = T)) |>
    add_row(data.frame(t = 0, prob = 0)) |>
    arrange(t) |>
    as.data.frame()
  
  intev        <- b - a
  idx          <- a <= indat_sum$t & indat_sum$t <= b
  time         <- indat_sum$t[idx]
  pred         <- indat_sum$prob[idx]
  # AUC using trapezoidal rule (pracma::trapz(time[ord],pred[ord])/intev)
  ord          <- order(time)
  area         <- sum(diff(time[ord])*zoo::rollmean(pred[ord], 2))/intev
  return(list(area = area, data = indat_sum)) 
}


#---------- Function to plot actual response
plotActualProb <- function(dat = ttrial, trial_var, resp_var, yaxis, lim) {
  
  indat <- as.data.frame(dat)
  indat$trial_var <- indat[, trial_var]
  indat$resp_var <- indat[, resp_var]
  
  plotactual <- indat |> filter(trial_var == 1) |>
    group_by(monthvisit_f, gr) |>
    summarise(prob = mean(resp_var, na.rm = T)) |>
    ggplot(aes(x = monthvisit_f, y = prob, group = gr, color = gr)) +
    geom_line(linewidth = 1) +
    scale_color_brewer(palette = "Set1", direction = -1) + 
    scale_y_continuous(limits = c(0, lim)) +
    labs(x = "Month", y = yaxis,
         color = NULL) +
    mytheme()
  
  return(plotactual)
}





#---------- Function to calculate area under the curve
auc_ttrial <- function(dat = ttrial, trial_var, a, b, outcome, boot = NULL, seed = NULL, conf_lv = 0.95) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  indat <- as.data.frame(dat)
  indat$outcome <- indat[, outcome]
  indat$trial_var <- indat[, trial_var]
  
  data_fit <- indat |> filter(trial_var == 1 & monthvisit > 0) |> 
    as.data.frame()
  
  #----- Parametric G-formula for analysis
  par_gform_ttrial <- function(dat, ...) {
    df <- dat
    gam1 <- gam(outcome ~ monthvisit_f + Tunder + ageunder +
                  sitenum + dep_treat + pde_yes + sex_yes + phy_yes + vit_yes +
                  s(pid, bs = "re"),
                family = binomial, method = "REML",
                data = df[df$tr == 1, ])
    
    # gam1 <- glmer(outcome ~ monthvisit_f + Tunder + ageunder + 
    #               sitenum + dep_treat + pde_yes + sex_yes + phy_yes + vit_yes + 
    #               (1 | pid),
    #             family = binomial, 
    #             data = df[df$tr == 1, ])
    
    gam0 <- gam(outcome ~ monthvisit_f + Tunder + ageunder +
                  sitenum + dep_treat + pde_yes + sex_yes + phy_yes + vit_yes +
                  s(pid, bs = "re"),
                family = binomial, method = "REML",
                data = df[df$tr == 0, ])
    
    # gam0 <- glmer(outcome ~ monthvisit_f + Tunder + ageunder + 
    #                 sitenum + dep_treat + pde_yes + sex_yes + phy_yes + vit_yes + 
    #                 (1 | pid),
    #               family = binomial, 
    #               data = df[df$tr == 0, ])
    
    df$pred1 <- predict.gam(gam1, newdata = df |> mutate(tr = 1), type = "response")
    df$pred0 <- predict.gam(gam0, newdata = df |> mutate(tr = 0), type = "response")
    adj_AUC1 <- auc_smooth_trap_ttrial(a = 0, b = 12, dat = df, resvar = "pred1", tvar = "monthvisit")
    adj_AUC0 <- auc_smooth_trap_ttrial(a = 0, b = 12, dat = df, resvar = "pred0", tvar = "monthvisit")
    AUC_diff <- adj_AUC1$area - adj_AUC0$area
    df_out   <- rbind(adj_AUC1$data, adj_AUC0$data)
    df_out$arm <- rep(c("Testosterone", "Placebo"), each = 5)
    return(list(AUC1 = adj_AUC1$area, AUC0 = adj_AUC0$area, AUC_diff = AUC_diff, data = df_out))
  }
  
  par_gform_out <- par_gform_ttrial(dat = data_fit)
  
  
  if (!is.null(boot)) {
    
    AUC1_b <- vector(length = boot)
    AUC0_b <- vector(length = boot)
    AUC_diff_b <- vector(length = boot)
    data_b <- NULL
    
    
    sam_id <- data_fit$pid |> unique()
    sam_n <- length(sam_id)
    
    #--- Bootstrap data
    for (i in 1:boot) {
      idx      <- sample(sam_id, sam_n, replace = T)
      freq_id  <- table(idx)
      df_boot  <- NULL
      
      for(j in 1:max(freq_id)) {
        # Loop over repeated id
        temp_df <- data_fit[data_fit$pid %in% names(freq_id[freq_id %in% c(j:max(freq_id))]), ]
        temp_df$boot_id <- paste0(temp_df$pid, "_", j)
        df_boot <- rbind(df_boot, temp_df)
      }
      
      par_gform_out_b <- par_gform_ttrial(dat = df_boot)
      
      AUC1_b[i] <- par_gform_out_b$AUC1
      AUC0_b[i] <- par_gform_out_b$AUC0
      AUC_diff_b[i] <- par_gform_out_b$AUC_diff
      tempdat <- par_gform_out_b$data
      tempdat$iter <- i
      data_b <- rbind(data_b, tempdat)
    }
    par_gform_out <- list(AUC1 = par_gform_out$AUC1, AUC0 = par_gform_out$AUC0, 
                          AUC_diff = par_gform_out$AUC_diff, data = par_gform_out$data,
                          AUC1_b = quantile(AUC1_b, prob = c((1 - conf_lv)/2, (1 + conf_lv)/2)), 
                          AUC0_b = quantile(AUC0_b, prob = c((1 - conf_lv)/2, (1 + conf_lv)/2)), 
                          AUC_diff_b = quantile(AUC_diff_b, prob = c((1 - conf_lv)/2, (1 + conf_lv)/2)), 
                          data_b = data_b)
  }
  
  return(par_gform_out)
}




#---------- Function to take the output of `auc_ttrial` to plot the response curve with shaded areas
plotAdjuProb <- function(dat, yaxis, lim) {
  indat <- dat
  plotadj <- indat |> 
    ggplot(aes(x = as.factor(t), y = prob, group = arm, color = arm)) +
    geom_area(data = subset(indat, arm == "Testosterone"), alpha = 0.1, fill = "#fb6a4a",
              show.legend = F) +
    geom_area(data = subset(indat, arm == "Placebo"), alpha = 0.3, fill = "#6baed6",
              show.legend = F) +
    geom_line(linewidth = 1) +
    scale_y_continuous(limits = c(0, lim)) +
    scale_color_brewer(palette = "Set1", direction = -1) + 
    labs(x = "Month", y = yaxis,
         color = NULL) +
    mytheme()
  
  return(plotadj)
}




#---------- Apply to Testosterone trial
#===============================================================================
#----- Trial 1: Sexual functional trial

# Verify result from the paper
# ttrial |> select(pid, monthvisit, Activity_bsl, Activity, chg_Activity, chg_Activity_cmr) |> View()

# m1 <- lmer(chg_Activity ~ monthvisit_f + tr + Activity_bsl + Tunder + 
#              ageunder + sitenum + dep_treat + pde_yes + phy_yes + vit_yes + (1 | pid),
#            REML = TRUE,
#            data = ttrial |> filter(sex_yes == 1 & monthvisit > 0))
# summary(m1)
# confint(m1)


#----- PDQ-Q4 score
ttrial <- ttrial |> mutate(chg_Activity_cmr = ifelse(chg_Activity >= 0.6, 1, 0))

# Estimate with bootstrap
pdqq4 <- auc_ttrial(dat = ttrial, trial_var = "sex_yes", a = 0, b = 12, 
                    outcome = "chg_Activity_cmr", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "sex_yes", resp_var = "chg_Activity_cmr", 
               yaxis = "Probability for change from baseline of \nPDQ-Q4 score >= 0.6", lim = 0.6)

png("Outfigs/ttrial_PDQQ4.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = pdqq4$data, yaxis = "Probability for change from baseline of \nPDQ-Q4 score \u2265 0.6", lim = 0.6)
dev.off()

c(pdqq4$AUC1, pdqq4$AUC1_b)
c(pdqq4$AUC0, pdqq4$AUC0_b)
c(pdqq4$AUC_diff, pdqq4$AUC_diff_b)


# Sensitivity analysis with different cutoffs
# cutoff <- seq(0.1, 4, by = 0.2)
# AUC_sen <- matrix(ncol = 10, nrow = length(cutoff))
# AUC_sen[, 1] <- cutoff
# 
# for (i in 1:length(cutoff)) {
#   ttrial1 <- ttrial |> mutate(chg_Activity_cmr = ifelse(chg_Activity >= cutoff[i], 1, 0))
#   pdqq4_i <- auc_ttrial(dat = ttrial1, trial_var = "sex_yes", a = 0, b = 12,
#                       outcome = "chg_Activity_cmr", boot = 200, seed = 12345)
#   AUC_sen[i, 2] <- pdqq4_i$AUC1
#   AUC_sen[i, 3] <- pdqq4_i$AUC1_b[[1]]
#   AUC_sen[i, 4] <- pdqq4_i$AUC1_b[[2]]
# 
#   AUC_sen[i, 5] <- pdqq4_i$AUC0
#   AUC_sen[i, 6] <- pdqq4_i$AUC0_b[[1]]
#   AUC_sen[i, 7] <- pdqq4_i$AUC0_b[[2]]
# 
#   AUC_sen[i, 8] <- pdqq4_i$AUC_diff
#   AUC_sen[i, 9] <- pdqq4_i$AUC_diff_b[[1]]
#   AUC_sen[i, 10] <- pdqq4_i$AUC_diff_b[[2]]
# }

#--- Running parallel
cutoff <- seq(0.1, 3.4, by = 0.1)

sens_ttrial_pdq <- function(cutoff) {
  AUC_sen <- matrix(ncol = 10, nrow = 1)
  AUC_sen[, 1] <- cutoff

  ttrial1 <- ttrial |> mutate(chg_Activity_cmr = ifelse(chg_Activity >= cutoff, 1, 0))
  pdqq4_i <- auc_ttrial(dat = ttrial1, trial_var = "sex_yes", a = 0, b = 12,
                        outcome = "chg_Activity_cmr", boot = 200, seed = 12345)
  AUC_sen[1, 2] <- pdqq4_i$AUC1
  AUC_sen[1, 3] <- pdqq4_i$AUC1_b[[1]]
  AUC_sen[1, 4] <- pdqq4_i$AUC1_b[[2]]

  AUC_sen[1, 5] <- pdqq4_i$AUC0
  AUC_sen[1, 6] <- pdqq4_i$AUC0_b[[1]]
  AUC_sen[1, 7] <- pdqq4_i$AUC0_b[[2]]

  AUC_sen[1, 8] <- pdqq4_i$AUC_diff
  AUC_sen[1, 9] <- pdqq4_i$AUC_diff_b[[1]]
  AUC_sen[1, 10] <- pdqq4_i$AUC_diff_b[[2]]

  return(AUC_sen)
}


start_time <- Sys.time()
r_pdq <-  plyr::llply(cutoff, function(i) sens_ttrial_pdq(cutoff = i), .parallel = TRUE)
AUC_sen_pdq <- do.call(rbind, r_pdq)
elapsed_time <- Sys.time()-start_time
print(elapsed_time)



png("Outfigs/ttrial_PDQQ4_sens.png", units="in", width = 8, height = 6, res = 300)
AUC_sen_pdq |> as.data.frame() |>
  ggplot(aes(x = V1, y = V8)) +
  geom_line(linewidth = 1, color = "#02818a") +
  geom_ribbon(aes(ymin = V9, ymax = V10), alpha = 0.3, fill = "gold") +
  scale_y_continuous(limits = c(0, 0.25)) +
  scale_x_continuous(limits = c(0, 3.5), breaks = seq(0, 3.5, by = 0.5)) +
  labs(x = "CMR cutoffs", y = "AUC (Testosterone - Placebo)",
       title = "AUC differences for different CMR cutoffs") +
  mytheme()
dev.off()





#----- DISF_M_II sexual desire score
ttrial <- ttrial |> mutate(chg_SEXDES_DER_cmr = ifelse(chg_SEXDES_DER >= 5, 1, 0))

# Estimate with bootstrap
disf <- auc_ttrial(dat = ttrial, trial_var = "sex_yes", a = 0, b = 12, 
                    outcome = "chg_SEXDES_DER_cmr", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "sex_yes", resp_var = "chg_SEXDES_DER_cmr", 
               yaxis = "Probability for change from baseline of \nDISF-M-II score \u2265 5", lim = 0.6)

png("Outfigs/ttrial_DISF_M_II.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = disf$data, yaxis = "Probability for change from baseline of \nDISF-M-II score \u2265 5", lim = 0.6)
dev.off()

c(disf$AUC1, disf$AUC1_b)
c(disf$AUC0, disf$AUC0_b)
c(disf$AUC_diff, disf$AUC_diff_b)


#--- Running parallel
cutoff_disf <- seq(1, 13, by = 0.5)

sens_ttrial_disf <- function(cutoff) {
  AUC_sen <- matrix(ncol = 10, nrow = 1)
  AUC_sen[, 1] <- cutoff
  
  ttrial1 <- ttrial |> mutate(chg_SEXDES_DER_cmr = ifelse(chg_SEXDES_DER >= cutoff, 1, 0))
  res_auc <- auc_ttrial(dat = ttrial1, trial_var = "sex_yes", a = 0, b = 12,
                        outcome = "chg_SEXDES_DER_cmr", boot = 200, seed = 12345)
  AUC_sen[1, 2] <- res_auc$AUC1
  AUC_sen[1, 3] <- res_auc$AUC1_b[[1]]
  AUC_sen[1, 4] <- res_auc$AUC1_b[[2]]
  
  AUC_sen[1, 5] <- res_auc$AUC0
  AUC_sen[1, 6] <- res_auc$AUC0_b[[1]]
  AUC_sen[1, 7] <- res_auc$AUC0_b[[2]]
  
  AUC_sen[1, 8] <- res_auc$AUC_diff
  AUC_sen[1, 9] <- res_auc$AUC_diff_b[[1]]
  AUC_sen[1, 10] <- res_auc$AUC_diff_b[[2]]
  
  return(AUC_sen)
}


start_time <- Sys.time()
r_disf <-  plyr::llply(cutoff_disf, function(i) sens_ttrial_disf(cutoff = i), .parallel = TRUE)
AUC_sen_disf<- do.call(rbind, r_disf)
elapsed_time <- Sys.time()-start_time
print(elapsed_time)



png("Outfigs/ttrial_DISF_M_II_sens.png", units="in", width = 8, height = 4, res = 300)
AUC_sen_disf |> as.data.frame() |>
  ggplot(aes(x = V1, y = V8)) +
  geom_line(linewidth = 1, color = "#02818a") +
  geom_ribbon(aes(ymin = V9, ymax = V10), alpha = 0.3, fill = "gold") +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 13), breaks = seq(0, 13, by = 2)) +
  labs(x = "CMR cutoffs", y = "AUC (Testosterone - Placebo)",
       title = "AUC differences for different CMR cutoffs") +
  mytheme()
dev.off()



#----- IIEF erectile function score

ttrial <- ttrial |> mutate(chg_IIEFEF_DER_cmr = ifelse(chg_IIEFEF_DER >= 5, 1, 0))

# Estimate with bootstrap
iief <- auc_ttrial(dat = ttrial, trial_var = "sex_yes", a = 0, b = 12, 
                   outcome = "chg_IIEFEF_DER_cmr", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "sex_yes", resp_var = "chg_IIEFEF_DER_cmr", 
               yaxis = "Probability for change from baseline of \nIIEF erectile function score \u2265 5", lim = 0.6)

png("Outfigs/ttrial_IIEF.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = iief$data, yaxis = "Probability for change from baseline of \nIIEF erectile function score \u2265 5", lim = 0.6)
dev.off()

c(iief$AUC1, iief$AUC1_b)
c(iief$AUC0, iief$AUC0_b)
c(iief$AUC_diff, iief$AUC_diff_b)




#--- Running parallel
cutoff_iief <- seq(1, 12, by = 0.5)

sens_ttrial_iief <- function(cutoff) {
  AUC_sen <- matrix(ncol = 10, nrow = 1)
  AUC_sen[, 1] <- cutoff
  
  ttrial1 <- ttrial |> mutate(chg_IIEFEF_DER_cmr = ifelse(chg_IIEFEF_DER >= cutoff, 1, 0))
  res_auc <- auc_ttrial(dat = ttrial1, trial_var = "sex_yes", a = 0, b = 12,
                        outcome = "chg_IIEFEF_DER_cmr", boot = 200, seed = 12345)
  AUC_sen[1, 2] <- res_auc$AUC1
  AUC_sen[1, 3] <- res_auc$AUC1_b[[1]]
  AUC_sen[1, 4] <- res_auc$AUC1_b[[2]]
  
  AUC_sen[1, 5] <- res_auc$AUC0
  AUC_sen[1, 6] <- res_auc$AUC0_b[[1]]
  AUC_sen[1, 7] <- res_auc$AUC0_b[[2]]
  
  AUC_sen[1, 8] <- res_auc$AUC_diff
  AUC_sen[1, 9] <- res_auc$AUC_diff_b[[1]]
  AUC_sen[1, 10] <- res_auc$AUC_diff_b[[2]]
  
  return(AUC_sen)
}


start_time <- Sys.time()
r_iief <-  plyr::llply(cutoff_iief, function(i) sens_ttrial_iief(cutoff = i), .parallel = TRUE)
AUC_sen_iief <- do.call(rbind, r_iief)
elapsed_time <- Sys.time()-start_time
print(elapsed_time)



png("Outfigs/ttrial_IIEF_sens.png", units="in", width = 8, height = 4, res = 300)
AUC_sen_iief |> as.data.frame() |>
  ggplot(aes(x = V1, y = V8)) +
  geom_line(linewidth = 1, color = "#02818a") +
  geom_ribbon(aes(ymin = V9, ymax = V10), alpha = 0.3, fill = "gold") +
  scale_y_continuous(limits = c(0, 0.3)) +
  scale_x_continuous(limits = c(0, 13), breaks = seq(0, 13, by = 2)) +
  labs(x = "CMR cutoffs", y = "AUC (Testosterone - Placebo)",
       title = "AUC differences for different CMR cutoffs") +
  mytheme()
dev.off()




#----- Trial 2: Physical Function Trial 
#===============================================================================

#----- Primary outcome: increase of ≥50 m in 6-min walk test 
# Check numbers in the original paper (same)
ttrial |> filter(phy_yes == 1 & monthvisit > 0) |>
  group_by(monthvisit, arm) |>
  summarize(n = sum(!is.na(WLK_Distance_bin), na.rm = T),
            yes = sum(WLK_Distance_bin == 1, na.rm = T),
            prob = mean(WLK_Distance_bin, na.rm = T))

ttrial <- ttrial |> mutate(WLK_Distance_bin_new = ifelse(monthvisit == 0, 0, WLK_Distance_bin))

# Estimate with bootstrap
wt6m <- auc_ttrial(dat = ttrial, trial_var = "phy_yes", a = 0, b = 12, 
                   outcome = "WLK_Distance_bin_new", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "phy_yes", resp_var = "WLK_Distance_bin_new", 
               yaxis = "Increase of \u226550 m in \n6-min walk test", lim = 0.6)

png("Outfigs/ttrial2_wt6m.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = wt6m$data, yaxis = "Increase of \u226550 m in \n6-min walk test", lim = 0.6)
dev.off()

c(wt6m$AUC1, wt6m$AUC1_b)
c(wt6m$AUC0, wt6m$AUC0_b)
c(wt6m$AUC_diff, wt6m$AUC_diff_b)


#----- Increase of ≥8 in PF-10 score
# Check numbers in the original paper (same)
ttrial |> filter(phy_yes == 1 & monthvisit > 0) |>
  group_by(monthvisit, arm) |>
  summarize(n = sum(!is.na(pf10_bin), na.rm = T),
            yes = sum(pf10_bin == 1, na.rm = T),
            prob = mean(pf10_bin, na.rm = T))

ttrial <- ttrial |> mutate(pf10_bin_new = ifelse(monthvisit == 0, 0, pf10_bin))


# Estimate with bootstrap
pf10 <- auc_ttrial(dat = ttrial, trial_var = "phy_yes", a = 0, b = 12, 
                   outcome = "pf10_bin_new", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "phy_yes", resp_var = "pf10_bin_new", 
               yaxis = "Increase of \u22658 in \nPF-10 score", lim = 0.6)

png("Outfigs/ttrial2_pf10.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = pf10$data, yaxis = "Increase of \u22658 in \nPF-10 score", lim = 0.6)
dev.off()

c(pf10$AUC1, pf10$AUC1_b)
c(pf10$AUC0, pf10$AUC0_b)
c(pf10$AUC_diff, pf10$AUC_diff_b)





# Estimate with bootstrap (all participants)
wt6m_all <- auc_ttrial(dat = ttrial, trial_var = "all_part", a = 0, b = 12, 
                   outcome = "WLK_Distance_bin_new", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "all_part", resp_var = "WLK_Distance_bin_new", 
               yaxis = "Increase of \u226550 m in \n6-min walk test", lim = 0.6)

png("Outfigs/ttrial2_wt6m_all.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = wt6m_all$data, yaxis = "Increase of \u226550 m in \n6-min walk test", lim = 0.6)
dev.off()

c(wt6m_all$AUC1, wt6m_all$AUC1_b)
c(wt6m_all$AUC0, wt6m_all$AUC0_b)
c(wt6m_all$AUC_diff, wt6m_all$AUC_diff_b)





# Estimate with bootstrap (all participants)
pf10_all <- auc_ttrial(dat = ttrial, trial_var = "all_part", a = 0, b = 12, 
                   outcome = "pf10_bin_new", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "all_part", resp_var = "pf10_bin_new", 
               yaxis = "Increase of \u22658 in \nPF-10 score", lim = 0.6)

png("Outfigs/ttrial2_pf10_all.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = pf10$data, yaxis = "Increase of \u22658 in \nPF-10 score", lim = 0.6)
dev.off()

c(pf10_all$AUC1, pf10_all$AUC1_b)
c(pf10_all$AUC0, pf10_all$AUC0_b)
c(pf10_all$AUC_diff, pf10_all$AUC_diff_b)







#----- Trial 3: Vitality Trial 
#===============================================================================

ttrial <- ttrial |> mutate(FACIT_bin_new = ifelse(monthvisit == 0, 0, FACIT_bin))

# Check numbers in the original paper (same)
ttrial |> filter(vit_yes == 1 & monthvisit > 0) |>
  group_by(monthvisit, arm) |>
  summarize(n = sum(!is.na(FACIT_bin), na.rm = T),
            yes = sum(FACIT_bin == 1, na.rm = T),
            prob = mean(FACIT_bin, na.rm = T))


# Estimate with bootstrap
facit <- auc_ttrial(dat = ttrial, trial_var = "vit_yes", a = 0, b = 12, 
                   outcome = "FACIT_bin_new", boot = 1000, seed = 12345)

plotActualProb(dat = ttrial, trial_var = "vit_yes", resp_var = "FACIT_bin_new", 
               yaxis = "Increase of \u22654 in \nFACIT–Fatigue score", lim = 1)

png("Outfigs/ttrial3_facit.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = facit$data, yaxis = "Increase of \u22654 in \nFACIT–Fatigue score", lim = 1)
dev.off()

c(facit$AUC1, facit$AUC1_b)
c(facit$AUC0, facit$AUC0_b)
c(facit$AUC_diff, facit$AUC_diff_b)




# Among all participants
# Check numbers in the original paper (same)
ttrial |> filter(all_part == 1 & monthvisit > 0) |>
  group_by(monthvisit, arm) |>
  summarize(n = sum(!is.na(FACIT_bin), na.rm = T),
            yes = sum(FACIT_bin == 1, na.rm = T),
            prob = mean(FACIT_bin, na.rm = T))


# Estimate with bootstrap
facit_all <- auc_ttrial(dat = ttrial, trial_var = "all_part", a = 0, b = 12, 
                    outcome = "FACIT_bin_new", boot = 1000, seed = 12345)

png("Outfigs/ttrial3_facit_all.png", units="in", width = 10, height = 6, res = 300)
plotAdjuProb(dat = facit_all$data, yaxis = "Increase of \u22654 in \nFACIT–Fatigue score", lim = 1)
dev.off()

c(facit_all$AUC1, facit_all$AUC1_b)
c(facit_all$AUC0, facit_all$AUC0_b)
c(facit_all$AUC_diff, facit_all$AUC_diff_b)





























