source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

source(here("code/programs/schonfeld_plot.R")) ## adapted plot code for Scales Schoenfeld residuals
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "PHchecks")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("eczema", "psoriasis")

exposure <- XX[1]
outcome <- YY[1]
ii <- 1

pdf(paste0(here("out/PHchecks"), "/ph_checks_severity.pdf"), width = 8, height = 10)
par(mfrow=c(4,3))
for(exposure in XX){
  ABBRVexp <- substr(exposure, 1, 3)
  if(ABBRVexp == "ecz"){ncat = 3}else{ncat = 2}
  for (outcome in YY){
    .dib(paste0(exposure,"~",outcome))
    # load df_model -----------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod4_severity_modeldata.rds"))
    
    # schoenfeld  -------------------------------------------------------------
    sch_resid <- cox.zph(cox_fit, transform = 'identity', terms = FALSE)
    mediator_est <- broom::tidy(cox_fit, conf.int = T, conf.level = 0.95, exp = T) %>% slice(1:ncat)
    
    ## severity model
    for(xx in 1:3){
      if(ncat==2 & xx == 3){
        plot(c(range(sch_resid$x), rev(range(sch_resid$x))),
             c(rep(mediator_est[1, "conf.low"], 2),rep(mediator_est[1, "conf.high"], 2)),
             col = 0, lty = 0, bty = "n",xaxt = "n", yaxt = "n", ylab = "", xlab ="")
      }else{
        plot_title <- str_remove(mediator_est[xx, "term"], "severity")
        # plot_schonfeld(sch_resid[xx], col = "darkgreen", df = 5, 
        #                thin_points = TRUE, thin_prop = 0.025,
        #                thin_col = ggplot2::alpha(1, 0.2),
        #                lwd = 1.5, resid = T,
        #                xlab = "Time (in days)", ylab = "")
        # mtext(expression(hat(beta)(t) ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
        # mtext(LETTERS[ii], side=3, line=2, col=1, cex=1, font=2, adj = 0)
        # if(xx==1){mtext(paste0(str_to_title(exposure), " ~ ", str_to_title(outcome)), side=3, line=2, col=1, cex=0.7, font=2, adj = 0.5)}
        # ii <- ii+1
        # 
        plot_schonfeld(sch_resid[xx], col = "darkgreen", df = 5,
                       lwd = 1.5, resid = F, xlab = "Time (in days)", hr = T, ylab = "")
        mtext(expression(e^{hat(beta)(t)} ~ "for" ~ exposed), side = 2, padj = -2, cex = 0.7)
        polygon(c(range(sch_resid$x), rev(range(sch_resid$x))),
                c(rep(mediator_est[xx, "conf.low"], 2),rep(mediator_est[xx, "conf.high"], 2)),
                col = ggplot2::alpha(4, 0.2), lty = 0)
        lines(range(sch_resid$x), rep(mediator_est[xx, "estimate"], 2), col = 4, lty = 4)
        mtext(LETTERS[ii], side=3, line=2, col=1, cex=1, font=2, adj = 0, padj = -0.5)
        mtext(paste0(str_to_title(exposure), " ~ ", str_to_title(outcome)), side=3, line=2, col=1, cex=0.7, font=2, adj = 1, padj = 0.5)
        mtext(plot_title, side=3, line=2, col=1, cex=0.7, font=2, adj = 1, padj = 1.75)
        ii <- ii+1
      }
    }
  }
}
dev.off()
