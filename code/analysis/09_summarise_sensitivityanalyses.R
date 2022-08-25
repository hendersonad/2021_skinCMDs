source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

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

plot_df <- plot_df_main %>% 
  mutate(analysis = "Main") %>% 
  bind_rows(mutate(plot_df_ethnicity, analysis = "Includes ethnicity (2006-2019)")) %>% 
  bind_rows(mutate(plot_df_noghosts_1yr, analysis = "Consult < 1yr before entry")) %>% 
  bind_rows(mutate(plot_df_noghosts_3yr, analysis = "Consult < 3yrs before entry")) 

plot_df$a[plot_df$analysis == "Main"] <- 1 
plot_df$a[plot_df$analysis != "Main"] <- 0.7 

plot_df$exposure[plot_df$exposure=="Atopic eczema"] <- "Eczema" ## bit of renaming if necessary

pdf(here::here("out/analysis/forest_plot15_sens_analyses_v2.pdf"), width = 11, height = 12)
pd <- position_dodge(width = 0.75)
plot_both <- ggplot(plot_df, aes(y = analysis, x = hr, xmin = ciL, xmax = ciU, group = outcome, colour = outcome, alpha = a, label = text_to_plot)) +
  geom_point(aes(x = hr), position = pd, size = 3, shape = 1) +
  geom_errorbar(aes(xmin  = ciL, xmax = ciU), position = pd, width = 0.25) +
  geom_vline(xintercept = 1, lty=2, col = alpha(1,0.4)) +  
  geom_text(aes(x = 0.999),
            alpha = 1,
            position = pd,
            size = 3.6,
            #colour = 1,
            show.legend = FALSE, 
            hjust = 1) +
  scale_x_log10(breaks=seq(0.5,2,0.1),limits=c(0.8,1.4)) +
  scale_y_discrete(limits=rev) +
  facet_grid(rows = vars(model), cols = vars(exposure)) +
  guides(colour = guide_legend("Outcome"), 
         alpha = "none") +
  labs(y = "Hazard ratio", x = "Model", caption = "(n) HR [95% CI]") +
  scale_alpha_identity() +
  theme_ali() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
print(plot_both)
dev.off()
