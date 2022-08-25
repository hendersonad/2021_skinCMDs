source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)

YY <- c("depression", "anxiety")
XX <- c("eczema", "psoriasis")

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
    
    # load original model
    mod_severity_mediated <-
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
    
    mod_severity_confounded <-
      coxph(
        Surv(t, out) ~ severity + carstairs + cal_period + comorbid + cci + strata(setid),
        data = df_model
      ) 
    saveRDS(
      mod_severity_confounded,
      file = paste0(
        datapath,
        "out/models_data/",
        ABBRVexp,
        "_",
        outcome,
        "_mod10_severity_confounded_modeldata.rds"
      )
    )
    
    tidy_results <- function(model){
      mod4tab <- broom::tidy(model, conf.int = T, exponentiate = T, conf.level = 0.95)
      mod4tab %>% 
        filter(str_detect(string = term, pattern = "^severity.")) %>% 
        mutate(Y = paste0(outcome), X = paste0(exposure)) %>% 
        select(X, Y, term, estimate, conf.low, conf.high)
    }
    mod4results_mediator <- tidy_results(mod_severity_mediated)
    mod4results_confound <- tidy_results(mod_severity_confounded)
    severity_results <- severity_results %>% 
      bind_rows(mutate(mod4results_mediator, model = "mediator"),
                mutate(mod4results_confound, model = "confounder"))
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

pdf(here::here("out/analysis/forest_plot2_severity_confounder_and_mediator.pdf"), width = 6, height = 4) 
p1 <- ggplot(severity_results, aes(x = severity, y = estimate, ymin = conf.low, ymax = conf.high, group = outcome, colour = outcome)) +
  geom_point(position = pd, size = 3, shape = 1) +
  geom_errorbar(position = pd, width = 0.25) +
  geom_hline(yintercept = 1, lty=2) +  
  #ylim(c(0,NA)) +
  scale_y_log10(breaks=seq(0.5,2,0.1),position="left",limits=c(ybase, yheight)) +
  scale_x_discrete(limits=rev) +
  facet_grid(rows = vars(model), cols = vars(exposure), drop = TRUE, space = "free", scales = "free") +
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

tabout <- severity_results %>% 
  arrange(model, exposure, outcome, severity) %>% 
  dplyr::select(-model) %>% 
  mutate(confint = paste0(round(conf.low, 2), " - ", round(conf.high, digits = 2))) %>% 
  dplyr::select(-conf.low, -conf.high) %>% 
  gt() %>% 
  tab_row_group(
    label = md("**Confounder adjusted model**"), 
    rows = 1:10
  ) %>% 
  tab_row_group(
    label = md("**Mediator adjusted model**"), 
    rows = 11:20
  ) %>% 
  fmt_number(columns = where(is.numeric), decimals = 2) %>% 
  cols_align(columns = 3:4, align = "right") %>% 
  cols_label(
    outcome = md("**Event**"),
    severity = md("**Eczema severity**"),
    estimate = md("**HR**"),
    confint = md("**95% CI**")
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = outcome)
  ) %>% 
  tab_footnote(
    footnote = "HR: Hazard ratio, CI: Confidence Interval",
    locations = cells_column_labels(
      columns = c(estimate, confint))) %>% 
  tab_spanner(
    label = md("**Severity model**"),
    columns = 3:4
  ) %>% 
  tab_footnote(
    footnote = "Additionally adjusted for BMI, alcohol misuse, smoking status and (eczema only) sleep problems and steroid use",
    locations = cells_column_spanners(
      spanners = "**Severity model**")) 

tabout %>%
  gt::gtsave(
    filename =  paste0("tab7b_regressionseverity_confounder_and_mediator.html"),
    path = here::here("out/tables")
  )
