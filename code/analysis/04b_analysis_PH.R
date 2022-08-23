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
ii <- 1

pdf(paste0(here("out/PHchecks"), "/ph_checks.pdf"), width = 8, height = 10)
par(mfrow=c(4,2))
for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  for (outcome in YY){
    .dib(paste0(exposure,"~",outcome))
    # load df_model -----------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    cox_fit3 <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod3_modeldata.rds"))
    
    # schoenfeld  -------------------------------------------------------------
    sch_resid3 <- cox.zph(cox_fit3, transform = 'identity')
    mediator_est <- broom::tidy(cox_fit3, conf.int = T, conf.level = 0.95, exp = T) %>% slice(1)
    
    ## mediator adjusted model
    plot_schonfeld(sch_resid3[1], col = "darkgreen", df = 5, 
                   thin_points = TRUE, thin_prop = 0.025,
                   thin_col = ggplot2::alpha(1, 0.2),
                   lwd = 1.5, resid = T,
                   xlab = "Time (in days)", ylab = "",
                   main = paste0(str_to_title(exposure), " ~ ", str_to_title(outcome)))
    mtext(expression(hat(beta)(t) ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    mtext(LETTERS[ii], side=3, line=2, col=1, cex=1, font=2, adj = 0)
    ii <- ii+1
    
    plot_schonfeld(sch_resid3[1], col = "darkgreen", df = 5,
                   lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "",
                   main = "Smoothed Schoenfeld residuals")
    mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
    polygon(c(range(sch_resid3$x), rev(range(sch_resid3$x))),
            c(rep(mediator_est$conf.low, 2),rep(mediator_est$conf.high, 2)),
            col = ggplot2::alpha(4, 0.2), lty = 0)
    lines(range(sch_resid3$x), rep(mediator_est$estimate, 2), col = 4, lty = 4)
    mtext(LETTERS[ii], side=3, line=2, col=1, cex=1, font=2, adj = 0)
    ii <- ii+1
  }
}
dev.off()