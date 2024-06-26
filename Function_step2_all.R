library(mgcv)
library(gratia)
library(numDeriv)
library(pROC)


# Set up data
#===============================================================================
df_wide_sim <- sim_dat_nocov(N = 100, n_measure = 25, diff_ratio = 0.08, 
                             arm_name = c("Treatment", "Control"))
df <- sim_dat_change_long_nocov(cutoff = -0.5, dat = df_wide_sim)

x <- plot_GRO(df)
x$step
x$line



auc_res <- auc_funct(a = 1, b = 25, dat = df)
auc_res2 <- auc_funct_idv(a = 1, b = 25, dat = df)
auc_res
auc_res2


#------ Smooth function using GAM
#===============================================================================
gam_m <- gam(prob ~ group + s(t, k = 15) + s(id, t, bs = "re"),
             data = df, method = "REML", family = binomial) 

summary(gam_m)
gam.check(gam_m)


# check actual vs prediction prob
df2 <- df
df2$pred <- predict(gam_m, type = "response")

df2 |> group_by(group, t) |>
  summarise(Observed = mean(prob),
            Predicted = mean(pred)) |>
  gather(-c(group, t), key = "type", value = prop) |>
  ggplot(aes(x = t, y = prop, color = group)) +
  geom_line(aes(linetype = type), linewidth = 1) +
  theme_bw()




#------ Function using trapezoidal rule
#===============================================================================

# Observed response prediction 
df$pred <- predict(gam_m, type = "response")
df$pred_var <- (predict(gam_m, type = "response", se.fit = TRUE)$se.fit)^2

# Counterfactual response prediction
df$pred_tr <- predict(gam_m, type = "response", newdata = df |> mutate(group = "Treatment"))
df$pred_ct <- predict(gam_m, type = "response", newdata = df |> mutate(group = "Control"))

# SE of prediction
df$pred_tr_var <- (predict(gam_m, type = "response", newdata = df |> mutate(group = "Treatment"), 
                           se.fit = TRUE)$se.fit)^2
df$pred_ct_var <- (predict(gam_m, type = "response", newdata = df |> mutate(group = "Control"), 
                           se.fit = TRUE)$se.fit)^2


auc_smooth_trap <- function(a, b,  dat, resvar, errvar, tvar = "t") {
  indat <- dat
  
  indat$pred <- indat[, resvar]
  indat$varfit <- indat[, errvar]
  indat$t <- indat[, tvar]
  
  indat_sum <- indat |> 
    group_by(t) |>
    summarise(prob = mean(pred),
              n = n(),
              var = (prob*(1-prob))/n + sum(varfit)/(n^2)) |>
    as.data.frame()
  
  intev        <- b - a
  idx          <- a <= indat_sum$t & indat_sum$t <= b
  time         <- indat_sum$t[idx]
  pred         <- indat_sum$prob[idx]
  var          <- indat_sum$var[idx]
  # AUC using trapezoidal rule (pracma::trapz(time[ord],pred[ord])/intev)
  ord          <- order(time)
  area         <- sum(diff(time[ord])*zoo::rollmean(pred[ord], 2))/intev
  var_area     <- sum((diff(time[ord])^2)*(zoo::rollsum(var[ord], 2)/4))/intev^2
  se_area      <- sqrt(var_area)
  return(list(area = area, se = se_area))
}

auc_smooth_trap(a = 0, b = 24, dat = df |> filter(group == "Treatment"),
               resvar = "pred", errvar = "pred_var")

auc_smooth_trap(a = 0, b = 24, dat = df |> filter(group == "Control"),
               resvar = "pred", errvar = "pred_var")


auc_smooth_trap(a = 0, b = 24, dat = df,
               resvar = "pred_tr", errvar = "pred_tr_var")

auc_smooth_trap(a = 0, b = 24,  dat = df,
               resvar = "pred_ct", errvar = "pred_ct_var")





#------ Test bootstrap
#===============================================================================
auc_bootA <- vector(length = 100)
auc_bootB <- vector(length = 100)

start_t <- Sys.time()

for (i in 1:100) {
  idx <- sample(1:200, 200, replace = T)
  bb <- table(idx)
  
  boot <- NULL
  for(zzz in 1:max(bb)) {
    # Loop over repeated id
    cc <- df[df$id %in% names(bb[bb %in% c(zzz:max(bb))]), ]
    cc$bid <- paste0(cc$id, zzz)
    boot <- rbind(boot, cc)
  }
  gam_boot <- gam(prob ~ group + s(t, k = 15) + s(id, t, bs = "re"),
                  data = boot, method = "REML", family = binomial)
  
  # Observed response prediction 
  boot$pred <- predict(gam_boot, type = "response")
  boot$pred_var <- (predict(gam_boot, type = "response", se.fit = TRUE)$se.fit)^2
  
  auc_outA <- auc_smooth_trap(a = 0, b = 24, model = gam_boot, dat = boot |> filter(group == "Treatment"),
                            resvar = "pred", errvar = "pred_var")
  auc_outB <- auc_smooth_trap(a = 0, b = 24, model = gam_boot, dat = boot |> filter(group == "Control"),
                            resvar = "pred", errvar = "pred_var")
  
  auc_bootA[i] <- auc_outA$area
  auc_bootB[i] <- auc_outB$area
}
Sys.time() - start_t



mean(auc_bootA)
sd(auc_bootA)

mean(auc_bootB)
sd(auc_bootB)










