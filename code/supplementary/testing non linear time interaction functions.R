
library(mfp)
library(visreg)

data(GBSG)
GBSG %>% head()
testmfp <- mfp(Surv(rfst, cens) ~ fp(age)+prm+esm+tumsize+menostat+tumgrad+strata(htreat), 
               family = cox, data = GBSG, select = 0.05)

testmfp_summ <- testmfp %>% summary()
testmfp_summ$call

refitmfp_original <- coxph(formula = Surv(rfst, cens) ~ strata(htreat) + prm + tumsize + 
                             tumgrad + I((age/100)^-1) + I((age/100)^-1 * log((age/100))), 
                           data = GBSG)

visreg(refitmfp_original, "age", ylab = "log(HR)")

## repeat for my data
small_df_ids <- sample(unique(df_model$setid), size = 1000)
small_df <- df_model %>% 
  filter(setid %in% small_df_ids)

testmfp <- mfp(Surv(t, out) ~ fp(exposed) + strata(setid), 
               family = cox, data = small_df, select = 0.05)
testmfp_summ <- testmfp %>% summary()
testmfp_summ$call 

refitmfp <- coxph(formula = Surv(t, out) ~ strata(setid) + I((t/1000)^1) +  I((t/1000)^2) + exposed, 
                  data = small_df)

refitmfp %>% summary()

visreg(refitmfp, "t", ylab = "log(HR)")
visreg(refitmfp, "exposed", ylab = "log(HR)")



# none of this is working particularly well. Aalen model?  ----------------
library(timereg)

small_df_ids <- sample(unique(df_model$setid), size = 1000)
small_df <- df_model %>% 
  filter(setid %in% small_df_ids)
small_df %>% head()
mod.aalen <- aalen(Surv(t, out) ~ exposed, data = df_model)
plot(mod.aalen)





df_model$exp <- as.numeric(df_model$exposed)-1

## Therneau and time dependent variables 
cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod1_modeldata.rds"))
cox_fit %>% summary() ## exp(beta) = 1.274

zp <- cox.zph(cox_fit, transform = function(time) log(time + 20))
plot(zp, resid = F)

zpI <- cox.zph(cox_fit, transform = "identity")
plot(zpI, resid = F)
plot_schonfeld(zpI[1], col = "darkgreen", df = 5,
               thin_points = TRUE, thin_prop = 0.025,
               thin_col = ggplot2::alpha(1, 0.2),
               lwd = 1.5, resid = T,
               xlab = "Time (in days)", ylab = "",
               main = "Minimally adjusted model")
plot_schonfeld(zpI[1], col = "darkgreen", df = 5,
               lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")

terry <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid),
               data = df_model,
               tt = function(x, t, ...) x * log(t+20))
terry
plot(zp[1], col = 12, lwd = 2, resid = F)
abline(h = coef(cox_fit)[1], lwd = 2, lty = 1, col = 9)
abline(coef(terry)[1:2], lwd = 2, lty = 1, col = 2)

## try with pspline
pspline <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid), 
                 data = df_model, 
                 tt = function(x, t, ...) x * pspline(t/365.25))

pspline$coefficients[1] %>% exp()
pspline$coefficients[-1]
basis_fn <- length(pspline$coefficients[!is.na(pspline$coefficients)])
output <- data.frame(fup = seq(min(df_model$exp + df_model$t) / 365.25,
                               max(df_model$exp + df_model$t) / 365.25,
                               0.01))

pspline_time <- pspline(output$fup)
output$HR <- pspline$coefficients[1] + (pspline_time %*% pspline$coefficients[-1])

## try with splines::bs
kk <- c(1, 2, 3)
spline <- coxph(Surv(t, out) ~ exp + tt(exp) + strata(setid), 
                data = df_model, 
                tt = function(x, t, ...) x * splines::bs(t/365.25, degree = 3, knots = kk))

spline %>% summary()

spline$coefficients
basis_fn <- length(spline$coefficients[!is.na(spline$coefficients)])-1
output <- data.frame(fup = seq(min(df_model$exp) + min(df_model$t) / 365.25,
                               max(df_model$exp + df_model$t) / 365.25,
                               0.01))
spline_time <- splines::bs(output$fup, degree = 3, knots = kk)


output$HRspline <- spline$coefficients[1] + (spline_time %*% spline$coefficients[-1])

par(mfrow = c(2,1))
plot_schonfeld(zpI[1], col = "darkgreen", df = 5, se = F,
               lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")
mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
abline(h = exp(coef(cox_fit)[1]), lwd = 2, lty = 2, col = 2)
abline(h = 1, lwd = 2, lty = 2, col = 9)

plot(range(output$fup), range(c(exp(output$HR),exp(output$HRspline))), type = "n", log = "y",
     xlab="Time (in years)", ylab="")
mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
lines(output$fup, exp(output$HR))
lines(output$fup, exp(output$HRspline), col = 4)
abline(h = exp(coef(cox_fit)[1]), lwd = 2, lty = 2, col = 2)
abline(h = 1, lwd = 2, lty = 2, col = 9)