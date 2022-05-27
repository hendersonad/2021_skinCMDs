library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)
library(readstata13)
library(muhaz)
require(splines)
require(gplots)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

source(here("code/programs/schonfeld_plot.R")) ## adapted plot code for Scales Schoenfeld residuals
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "PHchecks")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

#exposure <- XX[1]
#outcome <- YY[1]

for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  for (outcome in YY){
    .dib(paste0(exposure,"~",outcome))
    # load df_model -----------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    
    cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod1_modeldata.rds"))
    cox_fit3 <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod3_modeldata.rds"))
    
    cox_model_unstrat <- coxph(Surv(t, out) ~ exposed, data = df_model)
    cox_fit_unstrat <- survfit(cox_model_unstrat, newdata = data.frame(exposed = c("Unexposed", str_to_title(exposure))))
    
    km_fit <- survfit(Surv(t, out) ~ exposed, data = df_model)
    
    pdf(paste0(here("out/PHchecks"), "/ph_checks_", ABBRVexp, "_", substr(outcome,1,3), ".pdf"), width = 10, height = 8)
    par(mfrow=c(3,2))
    
    # compare KM and Cox fits -------------------------------------------------
    plot(cox_fit_unstrat, fun = "s", col = c(2,4), ylab="S(t|x)", xaxt = "n", xlab = "Years")
    axis(side = 1, at = seq(0, max(df_model$t), 365.25), labels = seq(0, max(df_model$t/365.25), 1))
    par(new=T)
    plot(km_fit, fun = "s", col = c(2,4), lty = 2, xaxt = "n", yaxt = "n", ylab = "", xlab = "",
         main = paste0(str_to_title(exposure), " ~ ", str_to_title(outcome)), font = 2)
    legend("bottomleft",c("Unexposed, Cox","Exposed, Cox",
                          "Unexposed, KM","Exposed, KM"),
           col=rep(c(2,4),2),lty=c(1,1,2,2))
    mtext("A", side=3, line=2, col=1, cex=1, font=2, adj = 0)
    
    # plot log(-log(s)) for X=1 and X=0 ---------------------------------------
    maxT <- max(df_model$t)
    maxlogT <- max(df_model$t) %>% log() %>% ceiling()
    xdisplay <- seq(0,maxlogT,1)
    xused <- exp(xdisplay)
    plot(km_fit,fun="cloglog",xlab="Time in days (log scale)",ylab="log(-log S(t))",xaxt = "n",xlim=c(0.25,maxT),col=c(2,4))
    axis(side = 1, at = xused, labels = xdisplay)
    legend("bottomright",c("Unexposed","Exposed"),col=c(2,4),lty=1)
    mtext("B", side=3, line=2, col=1, cex=1, font=2, adj = 0)
    
    # test interaction with t -------------------------------------------------
    cox_test <- coxph(Surv(t, out) ~ exposed + exposed*t, data = df_model)
    interaction_test_val <- broom::tidy(cox_test, exp = T, conf.int = T, conf.level = 0.99) %>% slice(3)
    gamma <- interaction_test_val$estimate %>% signif(digits = 3)
    gamma_lci <- interaction_test_val$conf.low %>% signif(digits = 3)
    gamma_uci <- interaction_test_val$conf.high %>% signif(digits = 3)
    text_print <- paste0("Estimate for interaction with time (to 3 sig. dig.): ", gamma, " (", gamma_lci, " - ", gamma_uci, ")")
    text(x = exp(-1), y = log(-log(min(km_fit$surv)))-0.5, text_print, pos = 4)
    
    # schoenfeld  -------------------------------------------------------------
    sch_resid1 <- cox.zph(cox_fit, transform = 'identity')
    sch_resid3 <- cox.zph(cox_fit3, transform = 'identity')
    minimal_est <- broom::tidy(cox_fit, conf.int = T, conf.level = 0.99, exp = T) %>% slice(1)
    mediator_est <- broom::tidy(cox_fit3, conf.int = T, conf.level = 0.99, exp = T) %>% slice(1)
    
    ## simple minimally adjusted model
    plot_schonfeld(sch_resid1[1], col = "darkgreen", df = 5,
                   thin_points = TRUE, thin_prop = 0.025,
                   thin_col = ggplot2::alpha(1, 0.2),
                   lwd = 1.5, resid = T,
                   xlab = "Time (in days)", ylab = "",
                   main = "Minimally adjusted model")
    mtext(expression(hat(beta)(t) ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    mtext("C", side=3, line=2, col=1, cex=1, font=2, adj = 0)
    
    plot_schonfeld(sch_resid1[1], col = "darkgreen", df = 5,
                   lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")
    mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    polygon(c(range(sch_resid1$x), rev(range(sch_resid1$x))),
         c(rep(minimal_est$conf.low, 2),rep(minimal_est$conf.high, 2)),
         col = ggplot2::alpha(4, 0.2), lty = 0)
    lines(range(sch_resid1$x), rep(minimal_est$estimate, 2), col = 4, lty = 4)
    mtext("D", side=3, line=2, col=1, cex=1, font=2, adj = 0)
    
    ## mediator adjusted model
    plot_schonfeld(sch_resid3[1], col = "darkgreen", df = 5, 
                   thin_points = TRUE, thin_prop = 0.025,
                   thin_col = ggplot2::alpha(1, 0.2),
                   lwd = 1.5, resid = T,
                   xlab = "Time (in days)", ylab = "",
                   main = "Mediator adjusted model")
    mtext(expression(hat(beta)(t) ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    mtext("E", side=3, line=2, col=1, cex=1, font=2, adj = 0)
    
    plot_schonfeld(sch_resid3[1], col = "darkgreen", df = 5,
                   lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")
    mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    polygon(c(range(sch_resid3$x), rev(range(sch_resid3$x))),
            c(rep(mediator_est$conf.low, 2),rep(mediator_est$conf.high, 2)),
            col = ggplot2::alpha(4, 0.2), lty = 0)
    lines(range(sch_resid3$x), rep(mediator_est$estimate, 2), col = 4, lty = 4)
    mtext("F", side=3, line=2, col=1, cex=1, font=2, adj = 0)
  
    dev.off()
  }
}
