library(cmprsk)
library(survminer)
library(tidyverse)

# cumulative incidence ----------------------------------------------------
if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

for(exposure in XX) {
  ABBRVexp <- str_sub(exposure, 1 , 3)
  # load data ---------------------------------------------------------------
  df_model_anx <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_anxiety.rds"
    ))
  df_model_dep <-
    readRDS(paste0(
      datapath,
      "out/df_model",
      ABBRVexp,
      "_depression.rds"
    ))
  for(outcome in YY) {
    out <- substr(outcome, 1, 3)
      
    df_plot <- get(paste0("df_model_", out))  
    
    df_plot$plot_t = df_plot$t/365.25
    
    fit <-  survfit(Surv(plot_t, out) ~ as.factor(exposed), data = df_plot)
    assign(paste0("km_fit", ABBRVexp, out), fit)
  }
}

plot_function <- function(exposure, outcome, plotfit){
  col1 <- as.numeric(substr(str_to_upper(exposure), 1, 1)=="E") + 2
  plot(
    plotfit,
    conf.int = T,
    lty = 1:length(unique(df_plot$exposed)),
    fun = function(x) {1 - x},
    bty = "n",
    col = c(1,col1),
    main = paste0(str_to_title(exposure), " ~ " , str_to_title(outcome)),
    xlab = "Years", 
    ylab = "Cumulative Incidence",
    ylim = c(0, 0.12)
  )
  legend(0, 0.11, legend=c("Unexposed", str_to_title(exposure)),
         col=c(1,col1), lty=1:2, cex=0.8, bty = "n")
}
pdf(here::here("out/analysis/fig_cuminc.pdf"), width = 8, 6)
par(mfrow = c(2,2))
  plot_function("Eczema", "Anxiety", km_fiteczanx)
  plot_function("Eczema", "Depression", km_fiteczdep)
  plot_function("Psoriasis", "Anxiety", km_fitpsoanx)
  plot_function("Psoriasis", "Depression", km_fitpsodep)
dev.off()

