source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")
dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "tables")), showWarnings = FALSE)

for(exposure in XX) {
  ABBRVexp <- substr(exposure, 1 , 3)
  mod1_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod1_modeldata.rds"
    ))
  mod2_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod2_modeldata.rds"
    ))
  mod3_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod3_modeldata.rds"
    ))
  
  mod1_dep <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_depression_mod1_modeldata.rds"
      )
    )
  mod2_dep <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_depression_mod2_modeldata.rds"
      )
    )
  mod3_dep <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_depression_mod3_modeldata.rds"
      )
    )
  
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
  
  for(outcome in YY){
    .dib(paste0(outcome, "~", exposure))
    shortout <- substr(outcome, 1, 3)
    df_model <- get(paste0("df_model_", shortout))
    
    mod1 <- get(paste0("mod1_", shortout))
    mod2 <- get(paste0("mod2_", shortout))
    mod3 <- get(paste0("mod3_", shortout))
    
    mod1_tab <- mod1 %>%
      tbl_regression(exp = TRUE) %>%
      modify_footnote(update = estimate ~ "Adjusted for matched set (age, sex, GP)")
    
    mod2_full <- mod2 %>% 
      tbl_regression(exp = T, conf.level = 0.95) %>%
      modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
    
    mod3_full <- mod3 %>% 
      tbl_regression(exp = T, conf.level = 0.95) %>%
      modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
    
    if (exposure == "eczema") {
      mod3_full <- mod3_full %>%
        modify_footnote(
          update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status, sleep problems, oral glucocorticoid use"
        )
    } else if (exposure == "psoriasis") {
      mod3_full <- mod3_full %>%
        modify_footnote(update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status")
    }
    
    tbls_full <- tbl_merge(
      tbls = list(mod1_tab, mod2_full, mod3_full),
      tab_spanner = c("**Crude**", "**Confounder**", "**Mediator**")
    )
    tbls_full %>%
      gtsummary::as_gt() %>%
      gt::gtsave(
        filename =  paste0("tab11_", ABBRVexp, "_", outcome, "_full.html"),
        path = here::here("out/tables")
      )
  }
}
