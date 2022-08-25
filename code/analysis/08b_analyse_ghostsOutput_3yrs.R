source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

source(here("code/programs/fn_getplotdata.R")) # will need this little function later and in other codes

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(here("out", "data")), showWarnings = FALSE)


YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

exposure <- XX[1]
make_regression_tab <- function(exposure){
  ABBRVexp <- str_sub(exposure, 1 , 3)
  # load models -------------------------------------------------------------
  mod1_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod1_modeldata_noghosts-3yrs.rds"
    ))
  mod2_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod2_modeldata_noghosts-3yrs.rds"
    ))
  mod3_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod3_modeldata_noghosts-3yrs.rds"
    ))
  
  # load models -------------------------------------------------------------
  mod1_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod1_modeldata_noghosts-3yrs.rds"
    ))
  mod2_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod2_modeldata_noghosts-3yrs.rds"
    ))
  mod3_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod3_modeldata_noghosts-3yrs.rds"
    ))
  
  # load data ---------------------------------------------------------------
  df_model_anx <-
    readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_noghosts_3yrs.rds"))
  
  df_model_dep <-
    readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_noghosts_3yrs.rds"))
  
  # Get n, n events and person-years by group ------------------------------
  get_data_info <- function(model, data) {
    vars <- attr(terms(model), "term.labels")
    vars <- vars[1:length(vars) - 1]
    df <- data %>%
      dplyr::select(setid, patid, gender, age, pracid, out, all_of(vars), dob, indexdate, enddate, tstart, tstop, t) 
    dim(df)
    
    # convert to data.table for speed of selecting the last row by group
    dt <- setDT(df)
    dt <- dt[ complete.cases(dt) ]
    dt[, firstObsNum := tstart[1L], by = list(setid, patid)]
    dt[, full_t := tstop-firstObsNum]
    dt <- dt[ dt[order(setid, patid, tstart), .I[c(.N)], by = list(setid, patid)]$V1 ]
    
    # convert back to tibble for summary code
    dtib <- tibble(dt)
    
    dfsum <- dtib %>%
      group_by(exposed) %>%
      summarise(n = n(),
                nevents = sum(out),
                pyars_mil = sum(full_t/(365.25)) / 1e6) %>%
      ungroup()
    
    dfsum
  }
  
  mod1_desc_anx <- get_data_info(model = mod1_anx, data = df_model_anx)
  mod2_desc_anx <- get_data_info(model = mod2_anx, data = df_model_anx)
  mod3_desc_anx <- get_data_info(model = mod3_anx, data = df_model_anx)
  
  mod1_desc_dep <- get_data_info(mod1_dep, df_model_dep)
  mod2_desc_dep <- get_data_info(mod2_dep, df_model_dep)
  mod3_desc_dep <- get_data_info(mod3_dep, df_model_dep)
  
  # Get CIs and p-values ----------------------------------------------------
  getCI <- function(model,
                    sigF = 3,
                    ci_level = 0.95) {
    model_sum <- summary(model, conf.int = ci_level)
    paste0(signif(model_sum$conf.int[1, 3], sigF),
           "-",
           signif(model_sum$conf.int[1, 4], sigF))
  }
  
  models_list <- list(mod1_anx,
                      mod2_anx,
                      mod3_anx,
                      mod1_dep,
                      mod2_dep,
                      mod3_dep)
  
  ci_values <- sapply(models_list, getCI)
  names(ci_values) <- c("mod1_anx","mod2_anx","mod3_anx","mod1_dep","mod2_dep","mod3_dep")
  
  ## function to put together N, n_event/p-yars, HR, CI
  merge_results <- function(model = "mod1"){
    model_desc_anx <- get(paste0(model, "_desc_anx"))
    model_desc_dep <- get(paste0(model, "_desc_dep"))
    mod_anx <- get(paste0(model, "_anx"))
    mod_dep <- get(paste0(model, "_dep"))
    n <- c(
      " ",
      prettyNum(model_desc_anx$n, big.mark = ","),
      " ",
      prettyNum(model_desc_dep$n, big.mark = ",")
    )
    #col4 Nevent/pyars 
    events <- c(
      " ",
      paste0(prettyNum(model_desc_anx$nevents, big.mark = ","), "/", prettyNum(model_desc_anx$pyars_mil, digits = 3, big.mark = ",")),
      " ",
      paste0(prettyNum(model_desc_dep$nevents, big.mark = ","), "/", prettyNum(model_desc_dep$pyars_mil, digits = 3, big.mark = ","))
    )
    #col5 HR 
    hr <-
      c(
        " ",
        "-",
        prettyNum(
          exp(mod_anx$coefficients[paste0("exposed", str_to_title(exposure))]),
          digits = 3,
          big.mark = ","
        ),
        " ",
        "-",
        prettyNum(
          exp(mod_dep$coefficients[paste0("exposed", str_to_title(exposure))]),
          digits = 3,
          big.mark = ","
        )
      )
    #col8 CI 
    ci <- c(" ", "-", getCI(mod_anx), " ", "-", getCI(mod_dep))
    
    modelNum <- substr(model, 4, 4)
    outtibble <- tibble(n, events, hr, ci)
    names(outtibble) <- paste0(names(outtibble), modelNum)
    return(outtibble)
  }
  
  # regression table ---------------------------------------------------------
  #col1 exposure
  char <-
    c("  ",
      "  Unexposed",
      paste0("  ", str_to_title(exposure)),
      "  ",
      "  Unexposed",
      paste0("  ", str_to_title(exposure)))
  
  #col2 outcome
  out <- c("Anxiety", rep(" ", 2), "Depression", rep(" ", 2))
  
  out_table <-
    bind_cols(
      characteristic = char,
      outcome = out,
      merge_results("mod1"),
      merge_results("mod2"),
      merge_results("mod3")
    )
  out_table
}

pso_table <- make_regression_tab(XX[1])
ecz_table <- make_regression_tab(XX[2])

tab3 <- bind_rows(ecz_table, pso_table)

tab3_out <- tab3 %>%
  gt() %>% 
  tab_row_group(
    label = md("**Atopic eczema**"), 
    rows = 1:6
  ) %>% 
  tab_row_group(
    label = md("**Psoriasis**"), 
    rows = 7:12
  ) %>% 
  row_group_order(c("**Atopic eczema**", "**Psoriasis**")) %>% 
  cols_align(columns = 3:dim(tab3)[2], align = "right") %>% 
  cols_label(
    characteristic = md("**Exposure**"),
    outcome = md("**Event**"),
    n1 = md("**N (total)**"),
    events1 = md("**No. events/person-years (mil)**"),
    hr1 = md("**HR**"),
    ci1 = md("**95% CI**"),
    n2 = md("**N**"),
    events2 = md("**No. events/person-years (mil)**"),
    hr2 = md("**HR**"),
    ci2 = md("**95% CI**"),
    n3 = md("**N**"),
    events3 = md("**No. events/person-years (mil)**"),
    hr3 = md("**HR**"),
    ci3 = md("**95% CI**")
  ) %>% 
  tab_spanner(
    label = md("**Crude model**"),
    columns = 3:6
  ) %>% 
  tab_spanner(
    label = md("**Confounder model**"),
    columns = 7:10
  ) %>% 
  tab_spanner(
    label = md("**Mediator model**"),
    columns = 11:14
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = outcome)
  ) %>% 
  tab_footnote(
    footnote = "HR: Hazard ratio, CI: Confidence Interval",
    locations = cells_column_labels(
      columns = c(hr1, hr2, hr3, ci1, ci2, ci3))) %>% 
  tab_footnote(
    footnote = "Adjusted for matched set (age, sex, GP)",
    locations = cells_column_spanners(
      spanners = "**Crude model**")) %>% 
  tab_footnote(
    footnote = "Additionally adjusted for calendar period and comorbidities",
    locations = cells_column_spanners(
      spanners = "**Confounder model**")) %>% 
  tab_footnote(
    footnote = "Additionally adjusted for BMI, alcohol misuse, smoking status and (eczema only) sleep problems and steroid use",
    locations = cells_column_spanners(
      spanners = "**Mediator model**")
  )

tab3_out
tab3_out %>%
  gt::gtsave(
    filename =  paste0("tab14_regression_noghosts_3yrs.html"),
    path = here::here("out/tables")
  )

# make the plot -----------------------------------------------------------
plot_df <- make_df_forest()
saveRDS(plot_df, here::here("out/data/df_forest_noghosts-3yrs.rds"))

pdf(here::here("out/analysis/forest_plot8_noghosts_3yrs.pdf"), width = 6, height = 4)
pd <- position_dodge(width = 0.3)
plot_new <- ggplot(plot_df, aes(x = model, y = hr, ymin = ciL, ymax = ciU, group = exposure, colour = exposure, alpha = a)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(0.9,1.35)) +
  scale_x_discrete(limits=rev) +
  facet_wrap(~outcome, ncol = 1) +
  coord_flip() +
  guides(colour = guide_legend("Exposure"), 
         alpha = "none") +
  labs(y = "Hazard ratio", x = "Model") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(plot_new)
dev.off()

# make table summarising the consult  drop_out  ---------------------------
ecz_cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-eczema-main-mhealth.dta"))
ecz_nonghosts <- haven::read_dta(paste0(datapath, "out/variables-ecz-consultations-yrbeforeindex-3yrs.dta"))

pso_cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-psoriasis-main-mhealth.dta"))
pso_nonghosts <- haven::read_dta(paste0(datapath, "out/variables-pso-consultations-yrbeforeindex-3yrs.dta"))

get_table <- function(ABBRVexp) {
  data <- get(paste0(ABBRVexp, "_cohort"))
  ghosts_var <- get(paste0(ABBRVexp, "_nonghosts"))
  df_noghosts <- data %>% 
    left_join(ghosts_var, by = "patid") %>% 
    mutate(pre_cons = replace_na(consyrbeforeindex, 0))
  table_noghosts <- df_noghosts %>% 
    count(exposed, pre_cons) %>% 
    group_by(exposed) %>% 
    mutate(prop = prop.table(n)*100 %>% signif(digits = 2)) %>% 
    pivot_wider(exposed, names_from = pre_cons, values_from = c(n, prop))
  table_noghosts
}
ecz_tab <- get_table("ecz") %>% 
  mutate(exposed = case_when(
    exposed == 0 ~ "Unexposed",
    exposed == 1 ~ "Eczema"
    )
  )
pso_tab <- get_table("pso") %>% 
  mutate(exposed = case_when(
    exposed == 0 ~ "Unexposed",
    exposed == 1 ~ "Psoriasis"
    )
  )

ecz_tab$exposed <- factor(ecz_tab$exposed, levels = c("Unexposed", "Eczema"))
pso_tab$exposed <- factor(pso_tab$exposed, levels = c("Unexposed", "Psoriasis"))


tab_ghosts <- ecz_tab %>% 
  bind_rows(pso_tab)

gt_ghosts <- tab_ghosts %>%
  ungroup() %>% 
  dplyr::select(exposed, n_0, prop_0, n_1, prop_1) %>%
  gt::gt() %>%
  tab_row_group(label = "Eczema cohort",
                rows = 1:2) %>%
  tab_row_group(label = "Psoriasis cohort",
                rows = 3:4) %>%
  gt::tab_header(title = md("Proportion of people with/without a consultation < 3 yr before indexdate")) %>%
  gt::fmt_number(columns = c(3, 5), decimals = 1) %>%
  gt::fmt_number(columns = c(2, 4), decimals = 0) %>%
  gt::data_color(
    columns = c(prop_1),
    colors = scales::col_numeric(
      palette = paletteer::paletteer_c(palette = "viridis::inferno",
                                       n = 100) %>% as.character(),
      domain = c(0,100)
    )
  ) %>%
  cols_label(n_0 = md("*n* without a recent consultation"),
             n_1 = md("*n* with a recent consultation"),
             prop_0 = "%",
             prop_1 = "%") 
gt_ghosts

gt_ghosts %>% 
  gt::gtsave(
    filename =  paste0("ghosts_analysis_3yrs.html"),
    path = here::here("out/supplementary")
  )
