library(tidyverse)
library(data.table)
library(here)
library(magrittr)
library(survival)
library(gt)


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
dir.create(file.path(here("out", "tables")), showWarnings = FALSE)

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
  
  # load models -------------------------------------------------------------
  mod1_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod1_modeldata.rds"
    ))
  mod2_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod2_modeldata.rds"
    ))
  mod3_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod3_modeldata.rds"
    ))
  
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
  
  # Get n, n events and person-years by group ------------------------------
  get_data_info <- function(model, data) {
    vars <- attr(terms(model), "term.labels")
    vars <- vars[1:length(vars) - 1]
    df <- data %>%
      select(setid, patid, gender, age, pracid, out, all_of(vars), dob, indexdate, enddate, tstart, tstop, t) 
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
    filename =  paste0("tab6_regression.html"),
    path = here::here("out/tables")
  )

# make the plot -----------------------------------------------------------
get_plot_data <- function(pretty_table) {
  numeric_tab <- pretty_table %>%
    select(-starts_with("n")) %>%
    select(-starts_with("events")) %>%
    select(-characteristic,-outcome) %>%
    separate(ci1, c("ciL1", "ciU1"), sep = "-") %>%
    separate(ci2, c("ciL2", "ciU2"), sep = "-") %>%
    separate(ci3, c("ciL3", "ciU3"), sep = "-") %>%
    drop_na() %>% 
    filter(hr1 != "-") %>% 
    mutate_if(is.character, as.numeric) 
  
  plot_tab1 <- pretty_table %>%
    select(2) %>%
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
  mutate(model = factor(model, labels = c("Crude", "Confounder Adjusted", "Mediator adjusted")))

## add alpha parameter to grey out other models
plot_df$a <- 0.5
plot_df$a[plot_df$model == "Mediator adjusted"] <- 1

plot_df$ciU %>% max()

pd <- position_dodge(width = 0.3)
pdf(here::here("out/analysis/forest_plot1.pdf"), width = 6, height = 4)
p1 <- ggplot(plot_df, aes(x = model, y = hr, ymin = ciL, ymax = ciU, group = exposure, colour = exposure, alpha = a)) +
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
print(p1)
dev.off()
saveRDS(plot_df, here::here("out/data/df_forest_main.rds"))



# Estimates rate differences using HR outputs -----------------------------

XX <- c("psoriasis", "eczema")
exposure <- XX[1]
get_rate_difference <- function(exposure){
  
  ABBRVexp <- str_sub(exposure, 1 , 3)
  # load models -------------------------------------------------------------
  mod2_anx <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_anxiety_mod2_modeldata.rds"
    ))

  # load models -------------------------------------------------------------
  mod2_dep <-
    readRDS(paste0(
      datapath,
      "out/models_data/",
      ABBRVexp,
      "_depression_mod2_modeldata.rds"
    ))
  
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
  
  rd_calculation <- function(outcome){
    df_model <- get(paste0("df_model_", substr(outcome,1,3)))
    model_output <- get(paste0("mod2_", substr(outcome,1,3)))
    # get rid of patids with missing data
    df_model_noNA <- df_model %>% 
      dplyr::select(setid, patid, t, out, exposed, carstairs, cal_period, comorbid, cci) %>% 
      drop_na()
    
    ## get p-yars at risk
    ## get n events
    ## calculate incidence rate per 100,000 years
    exposed_data <- df_model_noNA %>% 
      filter(exposed==str_to_title(exposure))  
    exposed_summary <- exposed_data %>% 
      summarise(pyears = sum(t/365.25), events = sum(out)) %>% 
      mutate(pyears_000 = pyears/1e5,
             rate_exposed = events/pyears_000)
    exposed_rate <- exposed_summary$rate_exposed
    
    ## take confounder adjusted HR and invert it
    model_estimates <- model_output %>% 
      broom::tidy(conf.int = T) %>% 
      filter(str_detect("exposed", string = term)) %>% 
      mutate(HR = exp(estimate))
    invertedHR <- 1/model_estimates$HR
    
    ## estimate unexposed rate as exposed_rate * 1/HR
    unexposed_rate <- exposed_rate * invertedHR
    
    ## calc rate difference 
    rate_diff_simple <- exposed_rate - unexposed_rate
    
    ## bootstrap confidence interval 
    patids <- unique(exposed_data$patid) 
    n <- patids %>% length()
    b <- 5000
    
    exp_btsp <- vector(length = b)
    unexp_btsp <- vector(length = b)
    rd_btsp <- vector(length = b)
    for(btsp in 1:b){
      if(btsp%%1000 == 0){print(c("Run ", btsp))}
      hr_draw <- rnorm(1, mean = model_estimates$estimate,
                        sd = model_estimates$std.error) 
      index <- sample(patids, size = n, replace = TRUE)
      exposed_btsp <- exposed_data %>% 
        filter(patid %in% index)
      exposed_btsp_summ <- exposed_btsp %>% 
        summarise(pyears = sum(t/365.25), events = sum(out)) %>% 
        mutate(pyears_000 = pyears/1e5,
               rate_exposed = events/pyears_000)
      exp_btsp[btsp] <- exposed_btsp_summ$rate_exposed
      unexp_btsp[btsp] <- exp_btsp[btsp] * 1/exp(hr_draw)
      rate_diff[btsp] <- exp_btsp[btsp] - unexp_btsp[btsp]
    }
    
    sigD <- 1
    rate_exposed_CI <- paste(round(quantile(exp_btsp, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    rate_unexposed_CI <- paste(round(quantile(unexp_btsp, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    rate_diff_CI <- paste(round(quantile(rate_diff, prob = c(0.025, 0.975)),sigD), collapse = " - ")
    hr_CI <- paste(round(exp(c(model_estimates$conf.low, model_estimates$conf.high)), 2), collapse = " - ")
    
    data.frame(
      exp = exposure, 
      out = outcome,
      pyears = exposed_summary$pyears_000,
      events = exposed_summary$events, 
      rate_exposed = exposed_summary$rate_exposed,
      rate_exposed_CI = rate_exposed_CI,
      hr = model_estimates$HR,
      hr_CI = hr_CI,
      rate_unexposed = unexposed_rate,
      rate_unexposed_CI = rate_unexposed_CI,
      rate_diff = rate_diff_simple,
      rate_diff_CI = rate_diff_CI
    )
  }
  rd_anx <- rd_calculation("anxiety")
  rd_dep <- rd_calculation("depression")
  
  bind_rows(rd_anx, rd_dep)
}

pso_ratediff <- get_rate_difference("psoriasis")
ecz_ratediff <- get_rate_difference("eczema")

rate_differences <- pso_ratediff %>% 
  bind_rows(ecz_ratediff)

gt_rates <- rate_differences %>%
  dplyr::select(-exp) %>% 
  mutate(out = str_to_title(out)) %>% 
  gt::gt() %>%
  tab_row_group(label = "Eczema",
                rows = 3:4) %>%
  tab_row_group(label = "Psoriasis",
                rows = 1:2) %>%
  gt::fmt_number(columns = c(2,4,6,8,10), decimals = 1) %>%
  gt::fmt_number(columns = 3, decimals = 0) %>%
  cols_label(out = "Outcome",
             pyears = "Person-years (100,000)",
             events = "Events",
             rate_exposed = "Rate (exposed group)",
             hr = "Hazard ratio",
             rate_unexposed = "Rate (unexposed group)",
             rate_diff = "Rate difference"
             ) 

rate_differences %>% write_csv(here::here("out/supplementary/rate_difference.csv"))
gt_rates %>% gt::gtsave(filename = "rate_difference.html", path = here::here("out/supplementary/"))
