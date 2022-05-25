library(here)

source(here::here("code/analysis/05d_analyse_secondary_data.R"))
source(here::here("code/analysis/05e_full_regression_tables.R"))
source(here::here("code/analysis/07_analyse_ghosts_regressions.R"))
source(here::here("code/analysis/07b_analyse_ghostsOutput.R"))
source(here::here("code/analysis/08_analyse_ghosts_regressions_3yrs.R"))
source(here::here("code/analysis/08b_analyse_ghostsOutput_3yrs.R"))
source(here::here("code/analysis/09_analyse_imputed_regressions.R"))
source(here::here("code/analysis/09b_analyse_imputed_regressionsOutput.R"))
source(here::here("code/analysis/10_summarise_sensitivityanalyses.R"))
source(here::here("code/analysis/11_analysisPH_t_v2.R"))
source(here::here("code/analysis/12_predictive_models.R"))

# lapply(paste("package:", names(sessionInfo()$otherPkgs), sep=""), 
#        detach, 
#        character.only = TRUE, 
#        unload = TRUE)
