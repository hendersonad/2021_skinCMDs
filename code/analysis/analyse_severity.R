library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    #datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)


YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

severity_results <- NULL
st <- Sys.time()
for(exposure in XX){
  #exposure <- XX[1]
  ABBRVexp <- substr(exposure, 1 , 3)
  
  for (outcome in YY) {
    #outcome = YY[1]
    .dib(paste0(outcome, "~", exposure))
    
    # load models -------------------------------------------------------------
    mod4 <-
      readRDS(
        paste0(
          datapath,
          "out/models_data/",
          ABBRVexp,
          "_",
          outcome,
          "_mod4_severity_modeldata.rds"
        )
      )
    
    mod4tab <- broom::tidy(mod4, conf.int = T, exponentiate = T, conf.level = 0.99)
    mod4results <- mod4tab %>% 
      filter(str_detect(string = term, pattern = "^severity.")) %>% 
      mutate(Y = paste0(outcome), X = paste0(exposure)) %>% 
      select(X, Y, term, estimate, conf.low, conf.high)
    
    severity_results <- severity_results %>% 
      bind_rows(mod4results)
    
    #mod4_tab <- mod4 %>%
    #  tbl_regression(exp = T,
    #                 #include = "severity",
    #                 conf.level = 0.99)
    
  }
}
end <- Sys.time()
end-st

severity_results <- severity_results %>% 
  mutate_at("term", ~str_remove(., "^severity")) %>% 
  mutate_at(c("X", "Y"), ~str_to_title(.)) %>% 
  rename(severity = term, 
         exposure = X,
         outcome = Y) 


pd <- position_dodge(width = 0.3)

roundUp <- function(x) 10^ceiling(log10(x))
ybase <- -0.1 + severity_results$conf.low %>% min() %>% round(digits = 2) 
yheight <- 0.1 + severity_results$conf.high %>% max() %>% round(digits = 2) 

ggplot(severity_results, aes(x = severity, y = estimate, ymin = conf.low, ymax = conf.high, group = outcome, colour = outcome)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(ybase, yheight)) +
  scale_x_discrete(limits=rev) +
  facet_wrap(~exposure, ncol = 1) +
  coord_flip() +
  guides(colour = guide_legend("Outcome")) +
  labs(y = "Hazard ratio", x = "Exposure severity") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

dev.copy(pdf, here::here("out/analysis/forest_plot2_severity.pdf"), width = 6, height = 4); dev.off()