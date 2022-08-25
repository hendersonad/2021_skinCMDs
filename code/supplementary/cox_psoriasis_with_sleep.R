source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)

XX <- c("psoriasis")
YY <- c("anxiety", "depression")
samplingSmall <- F
export_plots <- F # will be written to F if samplingSmall == T
export_models <- T
export_tables <- F
small_sample_size <- 100

st_time <- Sys.time()
for (exposure in XX[1]) {
  ABBRVexp <- str_sub(exposure, 1 , 3)
  
  .dib(exposure)
 
  
  for (outcome in YY) {
    .dib(outcome)

    df_model <- readRDS(file = paste0(datapath, "out/models_data/df_modelpso_", outcome, ".rds"))    
    # Run a mediator Cox regression -------------------------------------------
    .dib("Running model 3 (mediator)")
    if (ABBRVexp == "ecz") {
      mod3 <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + gc90days + strata(setid),
          data = df_model
        ) 
    } else if (ABBRVexp == "pso") {
      mod3 <-
        coxph(
          Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + strata(setid),
          data = df_model
        ) 
    }
    
    if (export_models) {
      saveRDS(
        mod3,
        file = paste0(
          datapath,
          "out/models_data/",
          ABBRVexp,
          "_sleep_",
          outcome,
          "_mod3_modeldata.rds"
        )
      )
    }
  }
}
# anxiety -----------------------------------------------------------------
pso_anx <-  readRDS(
  file = paste0(
    datapath,
    "out/models_data/",
    "pso",
    "_sleep_",
    "anxiety",
    "_mod3_modeldata.rds"
  )
)
df_model <- readRDS(file = paste0(datapath, "out/models_data/df_modelpso_", "anxiety", ".rds"))    

mod3_full <- pso_anx %>% tbl_regression(exp = T, conf.level = 0.99) %>%
  modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
mod3_full

# depression -----------------------------------------------------------------


pso_dep <-  readRDS(
  file = paste0(
    datapath,
    "out/models_data/",
    "pso",
    "_sleep_",
    "depression",
    "_mod3_modeldata.rds"
  )
)
df_model <- readRDS(file = paste0(datapath, "out/models_data/df_modelpso_", "depression", ".rds"))    

mod3_full <- pso_dep %>% tbl_regression(exp = T, conf.level = 0.99) %>%
  modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
mod3_full
