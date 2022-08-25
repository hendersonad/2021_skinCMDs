source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))
exposure <- "eczema"
ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure)

cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
severity <- haven::read_dta(paste0(datapath, "out/variables-ecz-severity.dta"))

anxiety_cohort <- readRDS(paste0(datapath, "out/df_modelecz_anxiety.rds"))
depression_cohort <- readRDS(paste0(datapath, "out/df_modelecz_depression.rds"))

ecz_sleep_everything <- haven::read_dta(paste0(datapath, "out/variables-ecz-sleep-all-additionalinfo.dta"))
ecz_sleep_def <- haven::read_dta(paste0(datapath, "out/variables-ecz-sleep-definite.dta"))
ecz_sleep_all <- haven::read_dta(paste0(datapath, "out/variables-ecz-sleep-all.dta"))

# Get unique events in everything sleep -----------------------------------
glimpse(ecz_sleep_everything)
ecz_sleep_unique <- ecz_sleep_everything %>% 
  group_by(patid, eventdate, sysdate, prodcode, medcode) %>% 
  slice(1)

ecz_sleep_unique <- ecz_sleep_unique %>% 
  ungroup() %>% 
  dplyr::select(patid, eventdate, sysdate, src, Possibledrugs, prodcode, ingredient, medcode, readterm) %>% 
  mutate_at("Possibledrugs", ~ifelse(is.na(.), 0, .)) %>% 
  rename(poss = Possibledrugs)


## get baseline observations for both cohort
anxiety_baseline <- anxiety_cohort %>% 
  dplyr::select(setid, patid, indexdate, exposed, severity, sleep) %>% 
  group_by(setid, patid) %>% 
  slice(1)
depression_baseline <- depression_cohort %>% 
  dplyr::select(setid, patid, indexdate, exposed, severity, sleep) %>% 
  group_by(setid, patid) %>% 
  slice(1)

full_baseline <- anxiety_baseline %>% 
  bind_rows(depression_baseline)

full_baseline <- full_baseline %>% 
  distinct(setid, patid, .keep_all = TRUE)

dim(full_baseline)
dim(anxiety_baseline)
dim(depression_baseline)

plot_sleep_breakdowns <- function(cohort = a_dataset){
  baseline_sleep <- cohort %>% 
    left_join(ecz_sleep_unique, by = "patid")
  
  ## censor events that are _after_ the enddate for that patid
  baseline_sleep$censorevent <- as.numeric(baseline_sleep$eventdate>baseline_sleep$indexdate)
  baseline_sleep$censorevent[is.na(baseline_sleep$censorevent)] <- 0
  baseline_sleep$eventdate[baseline_sleep$censorevent == 1] <- NA
  
  ## censor eczema severity because it is weird and some unexposed end up with mod/severe eczema
  baseline_sleep$severity[baseline_sleep$exposed == "Unexposed"] <- "None"
  
  # run 2x2 tables ----------------------------------------------------------
  twoXtwo <- function(df, exp, out){
  	df1 <- df %>% 
  		ungroup() %>% 
  		dplyr::select(exp = {{ exp }}, out = {{ out }})
  	tab <- table(df1$exp, df1$out, useNA = "always")
  	tab_p <- prop.table(tab,1)
  	tibble(
  		exposure = exp,
  		val = rownames(tab),
  		No = tab[, 1],
  		No_pc = tab_p[, 1] * 100,
  		Yes = tab[, 2],
  		Yes_pc = tab_p[, 2] * 100,
  		Miss = tab[, 3]
  	)
  }
  
  df_out <- baseline_sleep
  
  ## will use src to classify as prescription or diagnosis
  # get rid of NA srcs
  df_out$src[is.na(df_out$src)] <- 0
  # replace src as NA if eventdate > indexdate
  df_out$src[df_out$censorevent==1] <- 0
  
  # repeat for poss 
  df_out$poss[is.na(df_out$poss)] <- 0
  df_out$poss[df_out$censorevent==1] <- 0
  
  df_tab <- df_out %>% 
    ungroup() %>% 
    mutate(readcode = ifelse(src == 1, 1, 0)) %>% 
    mutate(drug_poss = ifelse(src == 2, 1, 0)) %>% 
    mutate(drug_def = ifelse(src == 2 & poss == 0, 1, 0)) %>% 
    group_by(setid, patid) %>% 
    mutate(readcode = max(readcode),
           drug_poss = max(drug_poss),
           drug_def = max(drug_def)) %>% 
    slice(1)
  
  df_tab$sleep_def <- df_tab$readcode + df_tab$drug_def
  df_tab$sleep_def[df_tab$sleep_def == 2] <- 1
  
  # 2x2 tables for sleep  ---------------------------------------------------
  # by severity
  tab_sleep_readcode <- twoXtwo(df = df_tab, exp = "severity", out = "readcode") %>% drop_na() %>% mutate(sleep_class = "read")
  tab_sleep_drugs_all <- twoXtwo(df = df_tab, exp = "severity", out = "drug_poss") %>% drop_na() %>% mutate(sleep_class = "drugs_poss")
  tab_sleep_drugs_def <- twoXtwo(df = df_tab, exp = "severity", out = "drug_def") %>% drop_na() %>% mutate(sleep_class = "drugs_def")
  tab_sleep <- twoXtwo(df = df_tab, exp = "severity", out = "sleep_def") %>% drop_na() %>% mutate(sleep_class = "sleep_def")
  
  # overall exposed
  ecz_sleep_readcode <- twoXtwo(df = df_tab, exp = "exposed", out = "readcode") %>% drop_na() %>% mutate(sleep_class = "read")
  ecz_sleep_drugs_all <- twoXtwo(df = df_tab, exp = "exposed", out = "drug_poss") %>% drop_na() %>% mutate(sleep_class = "drugs_poss")
  ecz_sleep_drugs_def <- twoXtwo(df = df_tab, exp = "exposed", out = "drug_def") %>% drop_na() %>% mutate(sleep_class = "drugs_def")
  ecz_sleep <- twoXtwo(df = df_tab, exp = "exposed", out = "sleep_def") %>% drop_na() %>% mutate(sleep_class = "sleep_def")
  
  tab_sleep_out <- bind_rows(tab_sleep_readcode,
                             tab_sleep_drugs_def, 
                             tab_sleep_drugs_all,
                             tab_sleep,
                             ecz_sleep_readcode,
                             ecz_sleep_drugs_all,
                             ecz_sleep_drugs_def,
                             ecz_sleep) %>% 
    filter(val != "Unexposed")
  tab_sleep_out$val[tab_sleep_out$val=="Eczema"] <- "Eczema (any severity)"
  
  gt_sleep <- tab_sleep_out %>%
  	dplyr::select(-exposure,-Miss, -sleep_class) %>%
    mutate(val = str_to_title(val)) %>% 
  	gt::gt() %>%
    tab_row_group(label = "Definite sleep drugs + Read codes",
                  rows = c(20, 13:16)) %>%
    tab_row_group(label = "All sleep Read codes",
                  rows = c(17, 1:4)) %>%
    tab_row_group(label = "Definite sleep drugs",
                  rows = c(19, 5:8)) %>%
    tab_row_group(label = "All sleep drugs",
                  rows = c(18, 9:12)) %>%
    gt::tab_header(title = "Sleep problems by eczema severity") %>%
  	gt::fmt_number(columns = c(3, 5), decimals = 1) %>%
  	gt::fmt_number(columns = c(2, 4), decimals = 0) %>%
  	gt::data_color(
  		columns = c(Yes_pc),
  		colors = scales::col_numeric(
  			palette = paletteer::paletteer_c(palette = "viridis::inferno",
  																			 n = 100) %>% as.character(),
  			domain = c(
  				min(tab_sleep_out$Yes_pc, na.rm = T),
  				max(tab_sleep_out$Yes_pc, na.rm = T)
  			)
  		)
  	) %>%
    tab_footnote(footnote = md(paste0("*n* = ", prettyNum(dim(df_tab)[1], big.mark = ","))),
                locations = cells_title("title")) %>% 
  	cols_label(val = "Severity",
  						 No_pc = "%",
  						 Yes_pc = "%") 
  
  # gt_sleep %>% 
  # 	gt::gtsave(
  # 		filename =  paste0("eczema_sleep_cohort_v3.html"),
  # 		path = here::here("out/supplementary")
  # 	)
  
  # bubble plot for proportion with code ------------------------------------
  p3 <- tab_sleep_out %>% 
    filter(sleep_class != "drugs_poss") %>% 
    mutate(lab = paste0(prettyNum(Yes_pc, big.mark = ",", digits = 3), "%")) %>% 
    mutate(val = str_to_title(val)) %>% 
    mutate(val = factor(val, levels = c("None","Eczema (Any Severity)", "Mild", "Moderate", "Severe"))) %>% 
    mutate(sleep_class = 
             factor(sleep_class, 
                    levels = c("read", "drugs_poss", "drugs_def", "sleep_def"), 
                    labels = c("Diagnosis", "Prescriptions", "All prescriptions", "Diagnosis or prescription"))) %>%
    mutate(yes_mil = Yes/1e3) %>% 
  	ggplot(aes(x = val, y = sleep_class, size = yes_mil, fill = Yes_pc)) +
  	  geom_point(alpha = 0.6, shape = 21, color = 1) +
  	  geom_text(aes(label = lab), nudge_y = 0, size = 4, col = 1) +
  	  scale_size(range = c(10, 30), name = "n (thousand)") +
  	  scale_fill_viridis_c(option = "inferno", name = "*%* with code") +
  	  labs(x = "Eczema severity", y = "") +
  	  guides(size = guide_legend(title.theme = element_text(angle = 0, face = "italic"))) +
  	  scale_y_discrete(limits=rev) +
      guides(size = "none") + 
  	  theme_ali() +
  	  theme(legend.position = "top", 
  	        legend.title = element_markdown(),
  	        axis.line = element_line(size = rel(0.5)),
  	        panel.grid.major.y = element_blank(),
  	        panel.grid = element_blank())
  p3
  #dev.copy(pdf, here::here("out/supplementary/sleep_data.pdf"), width = 8, height = 6); dev.off()
  
  	
  # what ingredients and diagnoses are most prevalent  ----------------------
  df_codes <- df_out
  
  df_codes <- df_codes %>%
    ungroup() %>%
    mutate(readcode = ifelse(src == 1, readterm, NA)) %>%
    mutate(drug_poss = ifelse(src == 2, ingredient, NA)) %>%
    mutate(drug_def = ifelse(src == 2 & poss == 0, ingredient, NA)) 
  
  sleep_codes <- df_codes %>% 
    ungroup() %>% 
    dplyr::select(setid, patid, exposed, severity, readcode, drug_poss, drug_def) %>% 
    mutate(prescription_ingredient = drug_def) %>% 
    pivot_longer(cols = c(readcode, prescription_ingredient), names_to = "src") %>% 
    mutate(desc = str_remove(value, "\\[d\\]|c/o - |\\[x\\]"))
  
  ## replace any mention of insomnia (including nos, or possible) with "insomnia"
  sleep_codes$desc[str_detect(sleep_codes$desc, ".?insomnia.?")] <- "insomnia"
  
  ## get rid of duplicate diagnoses/ingredients in the same patid
  sleep_tab_unique <- sleep_codes %>% 
    ungroup() %>% 
    dplyr::select(setid, patid, exposed, severity, src, desc) %>% 
    distinct(setid, patid, desc, .keep_all = TRUE)
  
  sleep_tab <- sleep_tab_unique %>% 
    group_by(desc, src, severity) %>% 
    summarise(n = n()) %>% 
    filter(!is.na(desc)) %>% 
    arrange(-n) %>% 
    mutate(rank = fct_reorder(as.factor(desc), n)) 
  
  sleep_tab$rank <- fct_reorder(sleep_tab$rank, sleep_tab$n)
  
  sleep_tab$src <- factor(sleep_tab$src, 
                          levels = c("readcode", "prescription_ingredient"),
                          labels = c("Diagnosis", "Prescription"))
  
  # sleep_tab$exposed <- factor(sleep_tab$exposed,
  #                             labels = c("Unexposed", "With eczema"))
  # sleep_tab$bigN <- ifelse(sleep_tab$exposed=="Unexposed", n_un, n_exp)
  
  x <- data.frame(table(df_tab$severity))
  sleep_tab <- sleep_tab %>% left_join(x, by = c("severity" = "Var1"))
  
  sleep_tab$pc <- round((sleep_tab$n/sleep_tab$Freq)*100, 1)
  sleep_tab$prettyN <- paste0(prettyNum(sleep_tab$n, big.mark = ","))
  sleep_tab %>% 
    group_by(severity, src) %>% 
    summarise(sum(n))
  sleep_plot <- sleep_tab %>% 
    filter(pc > 0.1)
  
  p1 <- ggplot(sleep_plot, aes(x = rank, y = pc, group = src, fill = src)) +
    geom_col(alpha = 0.4, col = alpha(1, 0.5)) + 
    geom_text(aes(y = 10, label = prettyN), col = alpha(1, 0.7), hjust = "right") +
    coord_flip() +
    facet_wrap(~severity, ncol = 2) + 
    labs(x = "Sleep codes", 
         caption = "Codes in at least 0.1% of population",
         fill ="") + 
    scale_y_continuous(name="% with code at study entry", labels = scales::comma) + 
    theme_ali() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          legend.title = element_markdown(),
          axis.line = element_line(size = rel(0.5)),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(size = rel(0.5), color = rgb(0, 0, 0, 0.4)),
          panel.grid = element_blank())
  p1
  #dev.copy(pdf, here::here("out/supplementary/sleep_codes.pdf"), width = 8, height = 8); dev.off()
  
  sleep_tab_exp <- sleep_tab_unique %>% 
    group_by(desc, src, exposed) %>% 
    summarise(n = n()) %>% 
    filter(!is.na(desc)) %>% 
    arrange(-n) %>% 
    mutate(rank = fct_reorder(as.factor(desc), n)) 
  
  sleep_tab_exp$rank <- fct_reorder(sleep_tab_exp$rank, sleep_tab_exp$n)
  
  sleep_tab_exp$src <- factor(sleep_tab_exp$src, 
                          levels = c("readcode", "prescription_ingredient"),
                          labels = c("Diagnosis", "Prescription"))
  
  x <- data.frame(table(df_tab$exposed))
  sleep_tab_exp <- sleep_tab_exp %>% left_join(x, by = c("exposed" = "Var1"))
  
  sleep_tab_exp$exposed <- factor(sleep_tab_exp$exposed,
                             labels = c("Unexposed", "With eczema"))
  
  sleep_tab_exp$pc <- round((sleep_tab_exp$n/sleep_tab_exp$Freq)*100, 1)
  sleep_tab_exp$prettyN <- paste0(prettyNum(sleep_tab_exp$n, big.mark = ","), " (", sleep_tab_exp$pc, "%)")
  sleep_tab_exp %>% 
    group_by(exposed, src) %>% 
    summarise(sum(n))
  sleep_plot_exp <- sleep_tab_exp %>% 
    filter(pc > 0.1)
  
  p2 <- ggplot(sleep_plot_exp, aes(x = rank, y = pc, group = src, fill = src)) +
    geom_col(alpha = 0.4, col = alpha(1, 0.5)) + 
    geom_text(aes(y = 10, label = prettyN), col = alpha(1, 0.7), hjust = "right") +
    coord_flip() +
    facet_wrap(~exposed, ncol = 2) + 
    labs(x = "Sleep codes", 
         caption = "Codes in at least 0.1% of population",
         fill ="") + 
    scale_y_continuous(name="Prevalence", labels = scales::comma) + 
    theme_ali() +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          legend.title = element_markdown(),
          axis.line = element_line(size = rel(0.5)),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_line(size = rel(0.5), color = rgb(0, 0, 0, 0.4)),
          panel.grid = element_blank())
  p2
  #dev.copy(pdf, here::here("out/supplementary/sleep_codes_exposed.pdf"), width = 8, height = 8); dev.off()
  
  return(list(p1 = p1, p2 = p2, p3 = p3))
}

anxiety_plots <- plot_sleep_breakdowns(cohort = anxiety_baseline)
depression_plots <- plot_sleep_breakdowns(cohort = depression_baseline)

p_all_anx <- plot_grid(anxiety_plots$p3 + theme(plot.background = element_rect(color = 1)), anxiety_plots$p1 + theme(plot.background = element_rect(color = 1)), ncol = 1 , rel_heights = c(1, 1.5), labels = "AUTO")
p_exp_anx <- plot_grid(anxiety_plots$p3 + theme(plot.background = element_rect(color = 1)), anxiety_plots$p2 + theme(plot.background = element_rect(color = 1)), ncol = 1 , rel_heights = c(1, 1.5), labels = "AUTO")
p_all_dep <- plot_grid(depression_plots$p3 + theme(plot.background = element_rect(color = 1)), depression_plots$p1 + theme(plot.background = element_rect(color = 1)), ncol = 1 , rel_heights = c(1, 1.5), labels = "AUTO")
p_exp_dep <- plot_grid(depression_plots$p3 + theme(plot.background = element_rect(color = 1)), depression_plots$p2 + theme(plot.background = element_rect(color = 1)), ncol = 1 , rel_heights = c(1, 1.5), labels = "AUTO")

p1 <- anxiety_plots$p3 + labs(title = "Eczema ~ Anxiety cohort") + theme(plot.background = element_rect(color = 1), plot.title = element_text(face = 2, size = 16, hjust = 0.5))
p2 <- depression_plots$p3 + labs(title = "Eczema ~ Depression cohort") + theme(plot.background = element_rect(color = 1), plot.title = element_text(face = 2, size = 16, hjust = 0.5))
p3 <- anxiety_plots$p1 + theme(plot.background = element_rect(color = 1))
p4 <- depression_plots$p1 + theme(plot.background = element_rect(color = 1))

p_all <- plot_grid(p1, p2, p3, p4, 
                   nrow = 2, labels = "AUTO", 
                   rel_heights = c(1, 1.5), align = "v", axis = "r")
pdf(here::here("out", "supplementary", "sleep_codes_all.pdf"), width = 13, height = 10)
  p_all
dev.off()

pdf(here::here("out", "supplementary", "anxiety_sleep_codes_all.pdf"), width = 8, height = 12)
  p_all_anx
dev.off()
pdf(here::here("out", "supplementary", "anxiety_sleep_codes_exp.pdf"), width = 8, height = 12)
  p_exp_anx
dev.off()

pdf(here::here("out", "supplementary", "depression_sleep_codes_all.pdf"), width = 8, height = 12)
  p_all_dep
dev.off()
pdf(here::here("out", "supplementary", "depression_sleep_codes_exp.pdf"), width = 8, height = 12)
  p_exp_dep
dev.off()
