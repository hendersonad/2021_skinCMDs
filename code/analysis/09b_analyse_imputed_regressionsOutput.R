library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)
library(readstata13)

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
    readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_anxiety_noghosts-3yrs.rds"))
  
  df_model_dep <-
    readRDS(paste0(datapath, "out/models_data/df_model", ABBRVexp, "_depression_noghosts-3yrs.rds"))
  
  # Get n, n events and person-years by group ------------------------------
  get_data_info <- function(model, data) {
    vars <- attr(terms(model), "term.labels")
    vars <- vars[1:length(vars) - 1]
    df <- data %>%
      select(setid, patid, gender, age, pracid, out, all_of(vars), dob, indexdate = indexdate.x, enddate, tstart, tstop, t) %>%
      mutate(indexNum = as.numeric(indexdate-dob),
             full_t = tstop-indexNum) %>% 
      group_by(setid, patid) %>% 
      slice(n()) %>% 
      select(-tstart, -tstop, -t)
    
    dfsum <- df %>%
      group_by(exposed) %>%
      summarise(n = n(),
                nevents = sum(out),
                pyars_mil = sum(full_t/(365.25)) / 1e6) %>%
      ungroup()
    
    dfsum
  }
  
  mod1_desc_anx <- get_data_info(mod1_anx, df_model_anx)
  mod2_desc_anx <- get_data_info(mod2_anx, df_model_anx)
  mod3_desc_anx <- get_data_info(mod3_anx, df_model_anx)
  
  mod1_desc_dep <- get_data_info(mod1_dep, df_model_dep)
  mod2_desc_dep <- get_data_info(mod2_dep, df_model_dep)
  mod3_desc_dep <- get_data_info(mod3_dep, df_model_dep)
  
  # Get CIs and p-values ----------------------------------------------------
  getP <- function(model,
                   sigF = 3,
                   ci_level = 0.99) {
    model_sum <- summary(model, conf.int = ci_level)
    pval <-
      model_sum$coefficients[paste0("exposed", str_to_title(exposure)), 5] %>% signif(digits = 1)
    if (pval < 0.0001) {
      pvalout <- "***"
    } else if (pval < 0.001) {
      pvalout <- "**"
    } else if (pval < 0.01) {
      pvalout <- "*"
    } else {
      pvalout <- paste0(pval)
    }
    pvalout
  }
  getCI <- function(model,
                    sigF = 3,
                    ci_level = 0.99) {
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
  
  p_values <- sapply(models_list, getP)
  ci_values <- sapply(models_list, getCI)
  names(ci_values) <- names(p_values) <- c("mod1_anx",
                                           "mod2_anx",
                                           "mod3_anx",
                                           "mod1_dep",
                                           "mod2_dep",
                                           "mod3_dep")
  
  
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
  
  #col3 N
  colN <- c(
    " ",
    prettyNum(mod1_desc_anx$n, big.mark = ","),
    " ",
    prettyNum(mod1_desc_dep$n, big.mark = ",")
  )
  #col4 Nevent
  colevents <- c(
    " ",
    prettyNum(mod1_desc_anx$nevents, big.mark = ","),
    " ",
    prettyNum(mod1_desc_dep$nevents, big.mark = ",")
  )
  #col5 p-years
  colpyars <- c(
    " ",
    prettyNum(mod1_desc_anx$pyars_mil, digits = 1, big.mark = ",", ),
    " ",
    prettyNum(mod1_desc_dep$pyars_mil, digits = 1, big.mark = ",")
  )
  
  #col6 N (model 1)
  mod1N <-
    c(prettyNum(sum(mod1_desc_anx$n), big.mark = ","),
      rep(" ", 2),
      prettyNum(sum(mod1_desc_dep$n), big.mark = ","),
      rep(" ", 2))
  #col7 HR (model 1)
  mod1HR <-
    c(
      " ",
      "Ref",
      prettyNum(
        exp(mod1_anx$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      ),
      " ",
      "Ref",
      prettyNum(
        exp(mod1_dep$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      )
    )
  #col8 CI (model 1)
  mod1Ci <- c(" ", "-", ci_values[1], " ", "-", ci_values[4])
  #col9 p (model 1)
  mod1P <- c(" ", "-", p_values[1], " ", "-", p_values[4])
  
  #col10 N (model 2)
  mod2N <-
    c(prettyNum(sum(mod2_desc_anx$n), big.mark = ","),
      rep(" ", 2),
      prettyNum(sum(mod2_desc_dep$n), big.mark = ","),
      rep(" ", 2))
  #col11 HR (model 2)
  mod2HR <-
    c(
      " ",
      "Ref",
      prettyNum(
        exp(mod2_anx$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      ),
      " ",
      "Ref",
      prettyNum(
        exp(mod2_dep$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      )
    )
  #col12 CI (model 2)
  mod2Ci <- c(" ", "-", ci_values[2], " ", "-", ci_values[5])
  #col13 p (model 2)
  mod2P <- c(" ", "-", p_values[2], " ", "-", p_values[5])
  
  #col14 N (model 3)
  mod3N <-
    c(prettyNum(sum(mod3_desc_anx$n), big.mark = ","),
      rep(" ", 2),
      prettyNum(sum(mod3_desc_dep$n), big.mark = ","),
      rep(" ", 2))
  #col15 HR (model 3)
  mod3HR <-
    c(
      " ",
      "Ref",
      prettyNum(
        exp(mod3_anx$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      ),
      " ",
      "Ref",
      prettyNum(
        exp(mod3_dep$coefficients[paste0("exposed", str_to_title(exposure))]),
        digits = 3,
        big.mark = ","
      )
    )
  #col16 CI (model 3)
  mod3Ci <- c(" ", "-", ci_values[3], " ", "-", ci_values[6])
  #col17 p (model 3)
  mod3P <- c(" ", "-", p_values[3], " ", "-", p_values[6])
  
  out_table <-
    tibble(
      characteristic = char,
      outcome = out,
      n = colN,
      nevents = colevents,
      pyars = colpyars,
      n1 = mod1N,
      hr1 = mod1HR,
      ci1 = mod1Ci,
      p1 = mod1P,
      n2 = mod2N,
      hr2 = mod2HR,
      ci2 = mod2Ci,
      p2 = mod2P,
      n3 = mod3N,
      hr3 = mod3HR,
      ci3 = mod3Ci,
      p3 = mod3P
    )
  gt::gt(out_table) %>%
    gt::cols_align(columns = 3:dim(out_table)[2], align = "right")
}

pso_table <- make_regression_tab(XX[1])
ecz_table <- make_regression_tab(XX[2])

## combine table data
tab1 <- ecz_table$`_data`
tab2 <- pso_table$`_data`

tab3 <- bind_rows(tab1, tab2)

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
    n = md("**N (total)**"),
    nevents = md("**No. events**"),
    pyars = md("**Person-years (mil)**"),
    n1 = md("**N**"),
    hr1 = md("**HR**"),
    ci1 = md("**99% CI**"),
    p1 = md("***p***"),
    n2 = md("**N**"),
    hr2 = md("**HR**"),
    ci2 = md("**99% CI**"),
    p2 = md("***p***"),
    n3 = md("**N**"),
    hr3 = md("**HR**"),
    ci3 = md("**99% CI**"),
    p3 = md("***p***")
  ) %>% 
  tab_spanner(
    label = md("**Crude model**"),
    columns = 6:9
  ) %>% 
  tab_spanner(
    label = md("**Confounder model**"),
    columns = 10:13
  ) %>% 
  tab_spanner(
    label = md("**Mediator model**"),
    columns = 14:17
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
    footnote = "Additionally adjusted for BMI, alcohol misuse, smoking status",
    locations = cells_column_spanners(
      spanners = "**Mediator model**")) %>% 
  tab_footnote(
    footnote = "***: p<0.0001",
    locations = cells_column_labels(
      columns = c(p1,p2,p3))) 

tab3_out
tab3_out %>%
  gt::gtsave(
    filename =  paste0("table4_regression_noghosts-3yrs.html"),
    path = here::here("out/analysis")
  )
tab3_out %>%
  gt::gtsave(
    filename =  paste0("table4_regression_noghosts-3yrs.rtf"),
    path = here::here("out/analysis")
  )



# make the plot -----------------------------------------------------------
get_plot_data <- function(pretty_table) {
  numeric_tab <- pretty_table %>%
    select(-starts_with("p")) %>%
    select(-starts_with("n")) %>%
    select(-characteristic,-outcome) %>%
    separate(ci1, c("ciL1", "ciU1"), sep = "-") %>%
    separate(ci2, c("ciL2", "ciU2"), sep = "-") %>%
    separate(ci3, c("ciL3", "ciU3"), sep = "-") %>%
    mutate_if(is.character, as.numeric) %>%
    drop_na()
  
  plot_tab1 <- pretty_table %>%
    select(2) %>%
    filter(outcome != " ") %>%
    bind_cols(numeric_tab)
}

ecz_plot <- get_plot_data(ecz_table$`_data`) %>%
  mutate(exposure = "Atopic eczema")
pso_plot <- get_plot_data(pso_table$`_data`) %>%
  mutate(exposure = "Psoriasis")

plot_df <- bind_rows(ecz_plot, pso_plot)
plot_df <- plot_df %>% 
  pivot_longer(cols = starts_with(c("hr", "ci")),
               names_to = c("metric", "model"),
               names_pattern = "(.*)([0-9])") %>% 
  pivot_wider(id_cols = c(outcome, exposure, model), names_from = metric)

##  models as factors 
plot_df <- plot_df %>%
  mutate(model = factor(model, labels = c("Crude", "Confounder Adjusted", "Mediator adjusted")))

## add alpha parameter to grey out other models
plot_df$a <- 0.5
plot_df$a[plot_df$model == "Mediator adjusted"] <- 1

plot_df$ciU %>% max()

saveRDS(plot_df, here::here("out/data/df_forest_noghosts-3yrs.rds"))

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
dev.copy(pdf, here::here("out/analysis/forest_plot6_noghosts-3yrs.pdf"), width = 6, height = 4); dev.off()


# import the original main analysis forest plot data to merge ------------------
does_exist <- list.files(here::here("out/data"), "df_forest_main.rds") %>% length() %>% as.logical()

if(does_exist) {
  plot_df_main <- readRDS(here::here("out/data/df_forest_main.rds"))
  plot_df_noghosts <- plot_df
  
  plot_df <- plot_df_main %>% 
    mutate(analysis = "Main") %>% 
    bind_rows(mutate(plot_df_noghosts, analysis = "Consult < 3yrs before entry"))
  
  pd <- position_dodge(width = 0.3)
  plot_both <- ggplot(plot_df, aes(x = model, y = hr, ymin = ciL, ymax = ciU, group = outcome, colour = outcome, alpha = a)) +
    geom_point(position = pd, size = 3, shape = 1) +
    geom_errorbar(position = pd, width = 0.25) +
    geom_hline(yintercept = 1, lty=2) +  
    #ylim(c(0,NA)) +
    scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(0.9,1.35)) +
    scale_x_discrete(limits=rev) +
    facet_grid(rows = vars(analysis), cols = vars(exposure)) +
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
  dev.copy(pdf, here::here("out/analysis/forest_plot7_sens_mainVnoghosts-3yrs.pdf"), width = 6, height = 6); dev.off()
}


# make table summarising the consult  drop_out  ---------------------------
ecz_cohort <- readstata13::read.dta13(paste0(datapath, "out/getmatchedcohort-eczema-main-mhealth.dta"))
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

# Compare main/1-yr/3-yr --------------------------------------------------
# import the original main analysis forest plot data to merge ------------------
does_exist <- list.files(here::here("out/data"), "df_forest_noghosts.rds") %>% length() %>% as.logical()

if(does_exist) {
  plot_df_main <- readRDS(here::here("out/data/df_forest_main.rds"))
  plot_df_noghosts_1yr <- readRDS(here::here("out/data/df_forest_noghosts.rds"))
  plot_df_noghosts_3yr <- readRDS(here::here("out/data/df_forest_noghosts-3yrs.rds"))
  
  plot_df <- plot_df_main %>% 
    mutate(analysis = "Main") %>% 
    bind_rows(mutate(plot_df_noghosts_1yr, analysis = "Consult < 1yr before entry")) %>% 
    bind_rows(mutate(plot_df_noghosts_3yr, analysis = "Consult < 3yrs before entry")) 
  
  pd <- position_dodge(width = 0.3)
  plot_both <- ggplot(plot_df, aes(x = model, y = hr, ymin = ciL, ymax = ciU, group = outcome, colour = outcome, alpha = a)) +
    geom_point(position = pd, size = 3, shape = 1) +
    geom_errorbar(position = pd, width = 0.25) +
    geom_hline(yintercept = 1, lty=2) +  
    #ylim(c(0,NA)) +
    scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(0.9,1.35)) +
    scale_x_discrete(limits=rev) +
    facet_grid(rows = vars(analysis), cols = vars(exposure)) +
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
  dev.copy(pdf, here::here("out/analysis/forest_plot7_sens_mainVnoghosts-all.pdf"), width = 6, height = 6); dev.off()
}