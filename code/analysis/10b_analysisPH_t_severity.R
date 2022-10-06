source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

source(here("code/programs/schonfeld_plot.R")) ## adapted plot code for Scales Schoenfeld residuals
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "PHchecks")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

exposure <- XX[1]
outcome <- YY[1]

#pdf(paste0(here::here("out/analysis"), "/fig2_spline_time_estimates_withBootstrap.pdf"), 8, 8)
par(mfrow = c(2,2))
ii=1
for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1, 3)
  if(ABBRVexp == "ecz"){ncat = 3}else{ncat = 2}
  for(outcome in YY) {
    # load data ---------------------------------------------------------------
    df_model <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_", outcome,".rds"))
    .dib(exposure)
    
    ## Spline
    cox_fit <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod4_severity_modeldata.rds"))
    
    df_model$sev <- as.numeric(df_model$severity)-1 ## need a numeric exposure variable for tt() to work
    
    ## Therneau and time dependent variables 
    mediator_est <- broom::tidy(cox_fit, conf.int = T, conf.level = 0.95, exp = T) %>% slice(1:ncat)
    
    ## Bootstrap to get a confidence interval
    output <- data.frame(fup = seq(min(df_model$sev + df_model$t) / 365.25,
                                   max(df_model$sev + df_model$t) / 365.25,
                                   0.01))
    
    b <- 500
    data_ids <- unique(df_model$setid) 
    n <- length(data_ids)
    
    output_Pspline <- matrix(NA, nrow = b, ncol = length(output$fup))
    output_Bspline <- matrix(NA, nrow = b, ncol = length(output$fup))
    
    tstart <- Sys.time()
    set.seed(105351)
    for(btsp in 1:b){
      if(btsp %% 10 == 0){print(paste("run", btsp, "..."))}
      
      if(btsp == 1){ # for the first run just use the actual data 
        sample_data <- df_model
      }else{ # then sample and jumble up the data on each iteration of the bootstrap 
        sample <- sample(data_ids, size = n, replace = TRUE)
        sample_ids <- tibble(setid = sample) %>% arrange(setid)
        sample_data <- sample_ids %>% 
          left_join(df_model, by = "setid")
      }
      
      # bspline
      kk <- c(1, 2, 3)
        pspline_model <-
          coxph(
            Surv(t, out) ~ sev + tt(sev) + carstairs + cal_period + comorbid + cci + strata(setid),  ## won't work because treates severity as linear!
            data = sample_data,
            tt = function(x, t, ...) x * pspline(t/365.25)
          )
        bspline_model <-
          coxph(
            Surv(t, out) ~ sev + tt(sev) + carstairs + cal_period + comorbid + cci + strata(setid),
            data = sample_data,
            tt = function(x, t, ...) x * splines::bs(t/365.25, degree = 3, knots = kk)
          )
      pspline_time <- pspline(output$fup)
      pspline_time_coeffs <- pspline_model$coefficients[str_detect(string = names(pspline_model$coefficients), pattern = "ps\\(")]
      out_Pspline <- pspline_model$coefficients[1] + (pspline_time %*% pspline_time_coeffs)
      
      spline_time <- splines::bs(output$fup, degree = 3, knots = kk)
      spline_time_coeffs <- bspline_model$coefficients[str_detect(string = names(bspline_model$coefficients), pattern = "tt\\(exp")]
      out_Bspline <- bspline_model$coefficients[1] + (spline_time %*% spline_time_coeffs)
      
      output_Pspline[btsp, ] <- out_Pspline
      output_Bspline[btsp, ] <- out_Bspline
    }
    tstop <- Sys.time()
    print(tstop-tstart)
    
    medP <- apply(output_Pspline,2,function(x){median(x, na.rm=T)}) %>% exp()
    ciP1 <- apply(output_Pspline,2,function(x){quantile(x,0.025, na.rm=T)}) %>% exp()
    ciP2 <- apply(output_Pspline,2,function(x){quantile(x,0.975, na.rm=T)}) %>% exp()
    medB <- apply(output_Bspline,2,function(x){median(x, na.rm=T)}) %>% exp()
    ciB1 <- apply(output_Bspline,2,function(x){quantile(x,0.025, na.rm=T)}) %>% exp()
    ciB2 <- apply(output_Bspline,2,function(x){quantile(x,0.975, na.rm=T)}) %>% exp()
    
    main_Pspline <- output_Pspline[1,] %>% exp()
    main_Bspline <- output_Bspline[1,] %>% exp()
    
    ymax <- max(c(ciP2,ciB2))
    plot(range(output$fup), range(c(main_Pspline,main_Bspline)), type = "n",
         xlab="Time (in years)", ylab="", ylim = c(0.8, 1.4))
    lines(range(output$fup), rep(mediator_est$estimate,2), lwd = 2, lty = 1, col = 2)
    polygon(c(range(output$fup), rev(range(output$fup))),
            c(rep(mediator_est$conf.low, 2),rep(mediator_est$conf.high, 2)),
            col = ggplot2::alpha(2, 0.2), lty = 0)
    mtext(paste0("HR(t) for ", exposure), side = 2, padj = -4, cex = 0.9)
    # pspline
    #lines(output$fup, main_Pspline, col = 1)
    lines(output$fup, medP, col = 1, lty = 2)
    polygon(c(output$fup, rev(output$fup)),
            c(ciP1,rev(ciP2)),
            col = ggplot2::alpha(1, 0.2), lty = 0)
    #lines(output$fup, main_Bspline, col = 4)
    lines(output$fup, medB, col = 4, lty = 2)
    polygon(c(output$fup, rev(output$fup)),
            c(ciB1,rev(ciB2)),
            col = ggplot2::alpha(4, 0.2), lty = 0)
    
    abline(h = 1, lwd = 2, lty = 2, col = 9)
    legend(0, 0.99, legend = c("Proportional hazards", "Penalised spline", "Basis spline (1,2,3 years)"), col = c(2, 1,4), lty = 1, bty = "n")
    mtext(paste0(LETTERS[ii], ": ", str_to_title(exposure), " ~ ", str_to_title(outcome)), side=3, line=2, col=1, cex=1, font=2, adj = 0)
    ii <- ii+1
  }
}
dev.off()
# 
# pdf(paste0(here::here("out/analysis"), "/fig3_linear_time_estimates.pdf"), 8, 8)
# p1 <- cowplot::plot_grid(
#   interaction_ploteczanx,
#   interaction_ploteczdep,
#   interaction_plotpsoanx,
#   interaction_plotpsodep
# )
# print(p1)
# dev.off()
