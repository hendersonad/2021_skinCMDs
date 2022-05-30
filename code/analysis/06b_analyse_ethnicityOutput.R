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

exposure <- XX[1]
make_regression_tab <- function(exposure){
  ABBRVexp <- str_sub(exposure, 1 , 3)
  # load models -------------------------------------------------------------
  mod2_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod_ethnicity_confound.rds"
    ))
  mod3_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod_ethnicity_mod.rds"
    ))
  
  # load models -------------------------------------------------------------
  mod2_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod_ethnicity_confound.rds"
    ))
  mod3_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod_ethnicity_mod.rds"
    ))
  
  # load data ---------------------------------------------------------------
  df_model_anx <-
    readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety_2006on.rds"))
  
  df_model_dep <-
    readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression_2006on.rds"))
  
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
  
  mod2_desc_anx <- get_data_info(model = mod2_anx, data = df_model_anx)
  mod3_desc_anx <- get_data_info(model = mod3_anx, data = df_model_anx)
  
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
  
  models_list <- list(mod2_anx,
                      mod3_anx,
                      mod2_dep,
                      mod3_dep)
  
  ci_values <- sapply(models_list, getCI)
  names(ci_values) <- c("mod2_anx","mod3_anx","mod2_dep","mod3_dep")
  
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
    label = md("**Confounder model**"),
    columns = 3:6
  ) %>% 
  tab_spanner(
    label = md("**Mediator model**"),
    columns = 7:10
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = outcome)
  ) %>% 
  tab_footnote(
    footnote = "HR: Hazard ratio, CI: Confidence Interval",
    locations = cells_column_labels(
      columns = c(hr2, hr3, ci2, ci3))) %>% 
  tab_footnote(
    footnote = "Additionally adjusted for calendar period, comorbidities and ethnicity",
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
    filename =  paste0("tab12_regression_ethnicity.html"),
    path = here::here("out/tables")
  )

# make the plot -----------------------------------------------------------
get_plot_data <- function(pretty_table) {
  numeric_tab <- pretty_table %>%
    dplyr::select(-starts_with("n")) %>%
    dplyr::select(-starts_with("events")) %>%
    dplyr::select(-characteristic,-outcome) %>%
    separate(ci2, c("ciL2", "ciU2"), sep = "-") %>%
    separate(ci3, c("ciL3", "ciU3"), sep = "-") %>%
    drop_na() %>% 
    filter(hr2 != "-") %>% 
    mutate_if(is.character, as.numeric) 
  
  plot_tab1 <- pretty_table %>%
    dplyr::select(2) %>%
    filter(outcome != " ") %>%
    bind_cols(numeric_tab)
}

ecz_plot <- get_plot_data(ecz_table) %>%
  mutate(exposure = "Atopic eczema")
pso_plot <- get_plot_data(pso_table) %>%
  mutate(exposure = "Psoriasis")

plot_df <- bind_rows(ecz_plot, pso_plot)
plot_df <- plot_df %>% 
  pivot_longer(cols = starts_with(c("hr", "ci")),
               names_to = c("metric", "model"),
               names_pattern = "(.*)([0-9])") %>% 
  pivot_wider(id_cols = c(outcome, exposure, model), names_from = metric)

##  models as factors 
plot_df <- plot_df %>%
  mutate(model = factor(model, labels = c("Confounder Adjusted", "Mediator adjusted")))

## add alpha parameter to grey out other models
plot_df$a <- 0.5
plot_df$a[plot_df$model == "Mediator adjusted"] <- 1

plot_df$ciU %>% max()

saveRDS(plot_df, here::here("out/data/df_forest_ethnicity.rds"))

pdf(here::here("out/analysis/forest_plot13_ethnicity.pdf"), width = 6, height = 4)
pd <- position_dodge(width = 0.3)
plot_new <- ggplot(plot_df, aes(x = model, y = hr, ymin = ciL, ymax = ciU, group = exposure, colour = exposure, alpha = a)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(0.9,1.4)) +
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