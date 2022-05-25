library(tidyverse)
library(data.table)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)
library(readstata13)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(here("out", "data")), showWarnings = FALSE)


YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

# Compare main/1-yr/3-yr/imputed data --------------------------------------------------
# import the original main analysis forest plot data to merge ------------------
plot_df_main <- readRDS(here::here("out/data/df_forest_main.rds"))
plot_df_noghosts_1yr <- readRDS(here::here("out/data/df_forest_noghosts.rds"))
plot_df_noghosts_3yr <- readRDS(here::here("out/data/df_forest_noghosts-3yrs.rds"))
plot_df_ethnicity <- readRDS(here::here("out/data/df_forest_ethnicity.rds"))
plot_df_impute <- readRDS(here::here("out/data/df_forest_ethnicity.rds"))

plot_df <- plot_df_main %>% 
  mutate(analysis = "Main") %>% 
  bind_rows(mutate(plot_df_ethnicity, analysis = "Includes ethnicity (2006-2019)")) %>% 
  bind_rows(mutate(plot_df_noghosts_1yr, analysis = "Consult < 1yr before entry")) %>% 
  bind_rows(mutate(plot_df_noghosts_3yr, analysis = "Consult < 3yrs before entry")) %>% 
  bind_rows(mutate(plot_df_impute, analysis = "Impute missing data")) 

pd <- position_dodge(width = 0.3)
plot_both <- ggplot(plot_df, aes(x = analysis, y = hr, ymin = ciL, ymax = ciU, group = outcome, colour = outcome, alpha = a)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(0.9,1.4)) +
  scale_x_discrete(limits=rev) +
  facet_grid(rows = vars(model), cols = vars(exposure)) +
  coord_flip() +
  guides(colour = guide_legend("Exposure"), 
         alpha = "none") +
  labs(y = "Hazard ratio", x = "Model") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(plot_both)
dev.copy(pdf, here::here("out/analysis/forest_plot15_sens_analyses.pdf"), width = 6, height = 12); dev.off()
  