library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)

if (Sys.info()["user"] == "lsh1510922") {
  if (Sys.info()["sysname"] == "Darwin") {
    datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
  if (Sys.info()["sysname"] == "Windows") {
    datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
  }
}
source(here("code/programs/fn_getplotdata.R")) # will need this little function later and in other codes get_plot_data

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

    # load data ---------------------------------------------------------------
    df_model <-
      readRDS(paste0(
        datapath,
        "out/df_model",
        ABBRVexp,
        "_",
        outcome, 
        ".rds"
      ))
    
    
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
    
    mod4tab <- broom::tidy(mod4, conf.int = T, exponentiate = T, conf.level = 0.95)
    mod4results <- mod4tab %>% 
      filter(str_detect(string = term, pattern = "^severity.")) %>% 
      mutate(Y = paste0(outcome), X = paste0(exposure)) %>% 
      dplyr::select(X, Y, term, estimate, conf.low, conf.high)
    
    severity_results <- severity_results %>% 
      bind_rows(mod4results)
  }
}
  
severity_results <- severity_results %>% 
  mutate_at("term", ~str_remove(., "^severity")) %>% 
  mutate_at(c("X", "Y"), ~str_to_title(.)) %>% 
  rename(severity = term, 
         exposure = X,
         outcome = Y) 


pd <- position_dodge(width = 0.3)

ybase <- -0.1 + severity_results$conf.low %>% min() %>% round(digits = 2) 
yheight <- 0.1 + severity_results$conf.high %>% max() %>% round(digits = 2) 

pdf(here::here("out/analysis/forest_plot2_severity.pdf"), width = 6, height = 4) 
p1 <- ggplot(severity_results, aes(x = severity, y = estimate, ymin = conf.low, ymax = conf.high, group = outcome, colour = outcome)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(ybase, yheight)) +
  scale_x_discrete(limits=rev) +
  facet_grid(rows = vars(exposure), drop = TRUE, space = "free", scales = "free") +
  coord_flip() +
  guides(colour = guide_legend("Outcome")) +
  labs(y = "Hazard ratio", x = "Exposure severity") +
  scale_alpha_identity() +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
print(p1)
dev.off()


make_gt_results <- function(exposure){
  #exposure <- XX[1]
  ABBRVexp <- substr(exposure, 1 , 3)
  
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
  
  # load models -------------------------------------------------------------
  mod4_dep <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_depression_mod4_severity_modeldata.rds"
      )
    )
  mod4_anx <-
    readRDS(
      paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_anxiety_mod4_severity_modeldata.rds"
      )
    )
  
  # severity table with n,  events,  p-years --------------------------------
  get_data_info <- function(model, data) {
    vars <- attr(terms(model), "term.labels")
    vars <- vars[1:length(vars) - 1]
    df <- data %>%
      dplyr::select(setid, patid, gender, age, pracid, out, all_of(vars), t) %>%
      distinct(setid, patid, .keep_all = TRUE) %>%
      drop_na()
    dfsum <- df %>%
      group_by(severity) %>%
      summarise(n = n(),
                nevents = sum(out),
                pyars_mil = sum(t) / 1e6) %>%
      ungroup()
    # dfsum %>%
    # 	summarise_at(-1, ~sum(.)) %>%
    # 	bind_cols(exposed = "Total") %>%
    # 	bind_rows(dfsum)
    dfsum
  }
  
  mod4_desc_anx <- get_data_info(mod4_anx, df_model_anx)
  mod4_desc_dep <- get_data_info(mod4_dep, df_model_dep)
  
  # Get CIs and p-values ----------------------------------------------------
  getP <- function(model,
                   sigF = 3,
                   ci_level = 0.95) {
    model_sum <- summary(model, conf.int = ci_level)
    modelrownames <- rownames(model_sum$coefficients)
    sev_names <- modelrownames[str_detect(modelrownames, "severity")]
    pval <- model_sum$coefficients[sev_names, 5] %>% signif(digits = 1)
    pvalout <- vector()
    for(ii in 1:length(pval)){
      if (pval[ii] < 0.0001) {
        pvalout[ii] <- "***"
      } else if (pval[ii] < 0.001) {
        pvalout[ii] <- "**"
      } else if (pval[ii] < 0.01) {
        pvalout[ii] <- "*"
      } else {
        pvalout[ii] <- paste0(pval[ii])
      }
    }
    pvalout
  }
  getCI <- function(model,
                    sigF = 3,
                    ci_level = 0.95) {
    model_sum <- summary(model, conf.int = ci_level)
    modelrownames <- rownames(model_sum$coefficients)
    sev_names <- modelrownames[str_detect(modelrownames, "severity")]
    ciout <- vector()
    for(ii in 1:length(sev_names)){
      ciout[ii] <- paste0(signif(model_sum$conf.int[ii, 3], sigF),
             "-",
             signif(model_sum$conf.int[ii, 4], sigF))
    }
    ciout
  }
  
  models_list <- list(mod4_anx, mod4_dep)
  
  p_values <- sapply(models_list, getP)
  ci_values <- sapply(models_list, getCI)
  colnames(ci_values) <- colnames(p_values) <- c("mod4_anx", "mod4_dep")
  
  # regression table ---------------------------------------------------------
  sev_levels <- unique(df_model_anx$severity) %>% as.character()
  sev_levels <- sev_levels[sev_levels != "None"]
  n_levels <- length(sev_levels)
  
  #col1 severity
  char <- c("  ",
      "  Unexposed",
      paste0("  ", str_to_title(sev_levels)),
      "  ",
      "  Unexposed",
      paste0("  ", str_to_title(sev_levels)))
  
  #col2 outcome
  out <- c("Anxiety", rep(" ", n_levels+1), "Depression", rep(" ", n_levels+1))
  
  #col3 N
  colN <- c(
    " ",
    prettyNum(mod4_desc_anx$n, big.mark = ","),
    " ",
    prettyNum(mod4_desc_dep$n, big.mark = ",")
  )
  colevents <- c(
    " ",
    paste0(prettyNum(mod4_desc_anx$nevents, big.mark = ","), "/", prettyNum(mod4_desc_anx$pyars_mil, digits = 3, big.mark = ",")),
    " ",
    paste0(prettyNum(mod4_desc_dep$nevents, big.mark = ","), "/", prettyNum(mod4_desc_dep$pyars_mil, digits = 3, big.mark = ","))
  )
  #col7 HR (model 1)
  mod4HR <-
    c(
      " ",
      "Ref",
      prettyNum(
        exp(mod4_anx$coefficients[1:n_levels]),
        digits = 3,
        big.mark = ","
      ),
      " ",
      "Ref",
      prettyNum(
        exp(mod4_dep$coefficients[1:n_levels]),
        digits = 3,
        big.mark = ","
      )
    )
  #col8 CI (model 1)
  mod4Ci <- c(" ", "-", ci_values[,1], " ", "-", ci_values[,2])
  #col9 p (model 1)
  mod4P <- c(" ", "-", p_values[,1], " ", "-", p_values[,2])
  
  out_table <-
    tibble(
      characteristic = char,
      outcome = out,
      n = colN,
      events = colevents,
      hr4 = mod4HR,
      ci4 = mod4Ci,
      p4 = mod4P
    )
  gt::gt(out_table) %>%
    gt::cols_align(columns = 3:dim(out_table)[2], align = "right")
}

pso_table <- make_gt_results(XX[1])
ecz_table <- make_gt_results(XX[2])

tab1 <- ecz_table$`_data`
tab2 <- pso_table$`_data`

tab3 <- bind_rows(tab1, tab2)

tab3_out <- tab3 %>%
  gt() %>% 
  tab_row_group(
    label = md("**Atopic eczema**"), 
    rows = 1:10
  ) %>% 
  tab_row_group(
    label = md("**Psoriasis**"), 
    rows = 11:18
  ) %>% 
  row_group_order(c("**Atopic eczema**", "**Psoriasis**")) %>% 
  cols_align(columns = 3:dim(tab3)[2], align = "right") %>% 
  cols_label(
    characteristic = md("**Exposure severity**"),
    outcome = md("**Event**"),
    n = md("**N (total)**"),
    events = md("**No. events/person-years (mil)**"),
    hr4 = md("**HR**"),
    ci4 = md("**95% CI**"),
    p4 = md("***p***")
  ) %>% 
  tab_spanner(
    label = md("**Severity model**"),
    columns = 5:7
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = outcome)
  ) %>% 
  tab_footnote(
    footnote = "HR: Hazard ratio, CI: Confidence Interval",
    locations = cells_column_labels(
      columns = c(hr4, ci4))) %>% 
  tab_footnote(
    footnote = "Additionally adjusted for calendar period and comorbidities",
    locations = cells_column_spanners(
      spanners = "**Severity model**")) %>% 
  tab_footnote(
    footnote = "***: p<0.0001",
    locations = cells_column_labels(
      columns = c(p4))) 

tab3_out
tab3_out %>%
  gt::gtsave(
    filename =  paste0("tab7_regressionseverity.html"),
    path = here::here("out/tables")
  )

# make forest plot with HRs printed ---------------------------------------
get_n_ecz <- tab1 %>% 
  dplyr::select(characteristic, starts_with("n")) %>% 
  filter(n != " ") %>% 
  mutate_all(~str_remove_all(.,",")) %>% 
  mutate_at("n", as.numeric) %>% 
  mutate(outcome = c(rep("Anxiety",4), rep("Depression",4))) %>% 
  mutate_if(is.character, ~str_remove_all(., " ")) %>% 
  mutate(exposure = "Eczema")
get_n_pso <- tab2 %>% 
  dplyr::select(characteristic, starts_with("n")) %>% 
  filter(n != " ") %>% 
  mutate_all(~str_remove_all(.,",")) %>% 
  mutate_at("n", as.numeric) %>% 
  mutate(outcome = c(rep("Anxiety",3), rep("Depression",3))) %>% 
  mutate_if(is.character, ~str_remove_all(., " ")) %>% 
  mutate(exposure = "Psoriasis")

  
plot_df <- get_n_pso %>% 
  bind_rows(get_n_ecz) %>% 
  left_join(severity_results, by = c("outcome", "exposure", "characteristic" = "severity")) %>% 
  drop_na()
  
# convert results to string and pad so they have same width for printing on plots
plot_df$text_hr <- str_pad(round(plot_df$estimate,2), 4, pad = "0", side = "right")
plot_df$text_ciL <- str_pad(round(plot_df$conf.low,2), 4, pad = "0", side = "right")
plot_df$text_ciU <- str_pad(round(plot_df$conf.high,2), 4, pad = "0", side = "right")
plot_df$text_n <- prettyNum(plot_df$n, big.mark = ",")
plot_df$text_to_plot <- str_pad(paste0("(",
                                       plot_df$text_n,
                                       ") ",
                                       plot_df$text_hr, 
                                       " [", 
                                       plot_df$text_ciL,
                                       ",", 
                                       plot_df$text_ciU,
                                       "]"),
                                28, pad = " ", side = "left")

pdf(here::here("out/analysis/forest_plot2_severity_v2.pdf"), width = 7.5, height = 4.5)
pd <- position_dodge(width = 0.75)
p1 <- ggplot(plot_df, aes(y = characteristic, x = estimate, xmin = conf.low, xmax = conf.high, group = outcome, colour = outcome, label = text_to_plot)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_vline(xintercept = 1, lty=2, col = alpha(1,0.4)) +  
  geom_text(aes(x = 0.95),
            alpha = 1,
            position = pd,
            size = 3.6,
            #colour = 1,
            show.legend = FALSE, 
            hjust = 1) +
  scale_x_log10(breaks=seq(0.5,2,0.1),limits=c(0.8,NA)) +
  scale_y_discrete(limits=rev) +
  facet_grid(rows = vars(exposure), drop = TRUE, space = "free", scales = "free") +
  guides(colour = guide_legend("Outcome"), 
         alpha = "none") +
  labs(y = "Hazard ratio", x = "Exposure severity", caption = "(n) HR [95% CI]") +
  theme_ali() +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
print(p1)
dev.off()
