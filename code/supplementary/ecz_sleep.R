library(here)
library(flextable)
library(gt)
library(gtsummary)
library(janitor)
library(timetk)
library(skimr)
library(glue)
library(gridExtra)
library(tidyverse)
library(ggtext)
library(ggpubr)
library(cowplot)
library(patchwork)


if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

exposure <- "eczema"
ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure)

cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
severity <- haven::read_dta(paste0(datapath, "out/variables-ecz-severity.dta"))

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

cohort_severity <- cohort %>% 
  left_join(severity, by = "patid")
cohort_severity$modsevere2 <- as.numeric(cohort_severity$modsevere) + 1
cohort_severity$modsevere2[cohort_severity$exposed==1 & is.na(cohort_severity$modsevere)] <- 1
cohort_severity$modsevere2[is.na(cohort_severity$modsevere2)] <- 0

## get max severity value
cohort_severity <- cohort_severity %>% 
  group_by(setid, patid) %>% 
  mutate(severity = max(modsevere2, na.rm = T)) %>% 
  slice(n())

## convert to factor
cohort_severity$severity <- factor(cohort_severity$severity, 0:3, labels = c("none", "mild", "moderate", "severe"))
table(cohort_severity$severity)

cohort_sleep <- cohort_severity %>% 
  left_join(severity, by = "patid") %>% 
  dplyr::select(setid, patid, exposed, indexdate, enddate, severity) %>% 
  left_join(ecz_sleep_unique, by = "patid")

## censor events that are _after_ the enddate for that patid
cohort_sleep$censorevent <- as.numeric(cohort_sleep$eventdate>cohort_sleep$enddate)
cohort_sleep$censorevent[is.na(cohort_sleep$censorevent)] <- 0
cohort_sleep$eventdate[cohort_sleep$censorevent == 1] <- NA

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

df_out <- cohort_sleep

df_out$src[is.na(df_out$src)] <- 0
df_out$poss[is.na(df_out$poss)] <- 0

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
table(df_tab$severity)
table(cohort_severity$severity)

df_tab$sleep_def <- df_tab$readcode + df_tab$drug_def
df_tab[df_tab$sleep_def == 2] <- 1

# all sleep codes - drugs only --------------------------------------------
tab_sleep_readcode <- twoXtwo(df = df_tab, exp = "severity", out = "readcode") %>% drop_na() %>% mutate(sleep_class = "read")
tab_sleep_drugs_all <- twoXtwo(df = df_tab, exp = "severity", out = "drug_poss") %>% drop_na() %>% mutate(sleep_class = "drugs_poss")
tab_sleep_drugs_def <- twoXtwo(df = df_tab, exp = "severity", out = "drug_def") %>% drop_na() %>% mutate(sleep_class = "drugs_def")
tab_sleep <- twoXtwo(df = df_tab, exp = "severity", out = "sleep_def") %>% drop_na() %>% mutate(sleep_class = "sleep_def")

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
  filter(val != "0")
tab_sleep_out$val[tab_sleep_out$val=="1"] <- "Eczema (any severity)"

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

gt_sleep %>% 
	gt::gtsave(
		filename =  paste0("eczema_sleep_cohort_v3.html"),
		path = here::here("out/supplementary")
	)

# bubble plot for proportion with code ------------------------------------
p3 <- tab_sleep_out %>% 
	  mutate(lab = paste0(prettyNum(Yes_pc, big.mark = ",", digits = 3), "%")) %>% 
	  mutate(val = str_to_title(val)) %>% 
	  mutate(val = factor(val, levels = c("None", "Mild", "Moderate", "Severe", "Eczema (Any Severity)"))) %>% 
	  mutate(sleep_class = 
	           factor(sleep_class, 
	                  levels = c("read", "drugs_poss", "drugs_def", "sleep_def"), 
	                  labels = c("All sleep Read codes", "All sleep drugs", "Definite sleep drugs", "Definite sleep drugs + Read codes"))) %>% 
	ggplot(aes(x = val, y = sleep_class, size = Yes, fill = Yes_pc)) +
	  geom_point(alpha = 0.6, shape = 21, color = 1) +
	  geom_text(aes(label = lab), nudge_y = 0, size = 4, col = 1) +
	  scale_size(range = c(5, 20), name = "*n* with code") +
	  scale_fill_viridis_c(option = "inferno", name = "*%* with code") +
	  labs(x = "Eczema severity", y = "") +
	  guides(size = "none") +
	  scale_y_discrete(limits=rev) +
	  theme_ali() +
	  theme(legend.position = "top", 
	        legend.title = element_markdown(),
	        axis.line = element_line(size = rel(0.5)),
	        panel.grid.major.y = element_blank(),
	        panel.grid = element_blank())
p3
dev.copy(pdf, here::here("out/supplementary/sleep_data.pdf"), width = 8, height = 6); dev.off()

	
# what ingredients and diagnoses are most prevalent  ----------------------
df_codes <- cohort_sleep

df_codes$src[is.na(df_codes$src)] <- 0
df_codes$poss[is.na(df_codes$poss)] <- 0

# all read codes and prescription ingredients  --------------------------------------------
df_codetab <- df_codes %>% 
  ungroup() %>% 
  mutate(readcode = ifelse(src == 1, readterm, NA)) %>% 
  mutate(drug_poss = ifelse(src == 2, ingredient, NA)) %>% 
  mutate(drug_def = ifelse(src == 2 & poss == 0, ingredient, NA)) %>% 
  group_by(setid, patid) %>% 
  mutate(readcode = readcode[which(!is.na(readcode))[1]],
         drug_poss = drug_poss[which(!is.na(drug_poss))[1]],
         drug_def = drug_def[which(!is.na(drug_def))[1]]) %>% 
  slice(1)

sleep_codes <- df_codetab %>% 
  ungroup() %>% 
  dplyr::select(setid, patid, exposed, severity, readcode, drug_poss, drug_def) %>% 
  mutate(drug_combo = ifelse(is.na(drug_def), drug_poss, drug_def)) %>% 
  pivot_longer(cols = c(readcode, drug_combo), names_to = "src") %>% 
  mutate(desc = str_remove(value, "\\[d\\]|c/o - |\\[x\\]")) 

sleep_tab <- sleep_codes %>% 
  group_by(desc, src, severity) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(desc)) %>% 
  arrange(-n) %>% 
  mutate(rank = fct_reorder(as.factor(desc), n)) 

sleep_tab$rank <- fct_reorder(sleep_tab$rank, sleep_tab$n)

sleep_tab$src <- factor(sleep_tab$src, 
                        labels = c("Prescription", "Read code"))

# sleep_tab$exposed <- factor(sleep_tab$exposed,
#                             labels = c("Unexposed", "With eczema"))
# sleep_tab$bigN <- ifelse(sleep_tab$exposed=="Unexposed", n_un, n_exp)

x <- data.frame(table(cohort_severity$severity))
sleep_tab <- sleep_tab %>% left_join(x, by = c("severity" = "Var1"))

sleep_tab$pc <- round((sleep_tab$n/sleep_tab$Freq)*100, 1)
sleep_tab$prettyN <- paste0(prettyNum(sleep_tab$n, big.mark = ","))
sleep_tab %>% 
  group_by(severity, src) %>% 
  summarise(sum(n))
sleep_plot <- sleep_tab %>% 
  filter(pc > 0.5)

p1 <- ggplot(sleep_plot, aes(x = rank, y = pc, group = src, fill = src)) +
  geom_col(alpha = 0.4, col = alpha(1, 0.5)) + 
  geom_text(aes(y = 20, label = prettyN), col = alpha(1, 0.7), hjust = "right") +
  coord_flip() +
  facet_wrap(~severity, ncol = 2) + 
  labs(x = "Sleep codes", 
       caption = "Codes in at least 0.5% of population",
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
p1
dev.copy(pdf, here::here("out/supplementary/sleep_codes.pdf"), width = 8, height = 8); dev.off()



# sleep_tab$bigN <- ifelse(sleep_tab$exposed=="Unexposed", n_un, n_exp)

sleep_tab_exp <- sleep_codes %>% 
  group_by(desc, src, exposed) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(desc)) %>% 
  arrange(-n) %>% 
  mutate(rank = fct_reorder(as.factor(desc), n)) 

sleep_tab_exp$rank <- fct_reorder(sleep_tab_exp$rank, sleep_tab_exp$n)

sleep_tab_exp$src <- factor(sleep_tab_exp$src, 
                        labels = c("Prescription", "Read code"))

x <- data.frame(table(cohort_severity$exposed))
x$exposed <- as.numeric(x$Var1) - 1
sleep_tab_exp <- sleep_tab_exp %>% left_join(x, by = "exposed")

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
  geom_text(aes(y = 20, label = prettyN), col = alpha(1, 0.7), hjust = "right") +
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
dev.copy(pdf, here::here("out/supplementary/sleep_codes_exposed.pdf"), width = 8, height = 8); dev.off()

pdf(here::here("out/supplementary/sleep_codes_all.pdf"), width = 7, height = 10)
  plot_grid(p3, p1, ncol = 1 , rel_heights = c(1,2), labels = "AUTO")
dev.off()

pdf(here::here("out/supplementary/sleep_codes_exp.pdf"), width = 7, height = 10)
  plot_grid(p3, p2, ncol = 1 , rel_heights = c(1,2), labels = "AUTO")
dev.off()

# Eczema yes/no results ---------------------------------------------------
tab_sleep_readcode <- twoXtwo(df = df_tab, exp = "exposed", out = "readcode") %>% drop_na() %>% mutate(sleep_class = "read")
tab_sleep_drugs_all <- twoXtwo(df = df_tab, exp = "exposed", out = "drug_poss") %>% drop_na() %>% mutate(sleep_class = "drugs_poss")
tab_sleep_drugs_def <- twoXtwo(df = df_tab, exp = "exposed", out = "drug_def") %>% drop_na() %>% mutate(sleep_class = "drugs_def")
tab_sleep <- twoXtwo(df = df_tab, exp = "exposed", out = "sleep_def") %>% drop_na() %>% mutate(sleep_class = "sleep_def")

tab_sleep_out <- bind_rows(tab_sleep_readcode,
                           tab_sleep_drugs_def, 
                           tab_sleep_drugs_all,
                           tab_sleep)

gt_sleep <- tab_sleep_out %>%
  dplyr::select(-exposure,-Miss, -sleep_class) %>%
  mutate(val = str_to_title(val)) %>% 
  gt::gt() %>%
  tab_row_group(label = "Definite sleep drugs + Read codes",
                rows = 7:8) %>%
  tab_row_group(label = "All sleep Read codes",
                rows = 1:2) %>%
  tab_row_group(label = "Definite sleep drugs",
                rows = 3:4) %>%
  tab_row_group(label = "All sleep drugs",
                rows = 5:6) %>%
  gt::tab_header(title = "Sleep problems by eczema diagnosis") %>%
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
  cols_label(val = "Eczema",
             No_pc = "%",
             Yes_pc = "%") 

gt_sleep %>% 
  gt::gtsave(
    filename =  paste0("eczema_sleep_cohort_yesno.html"),
    path = here::here("out/supplementary")
  )

# bubble plot for proportion with code ------------------------------------
p3 <- tab_sleep_out %>% 
  mutate(lab = paste0(prettyNum(Yes_pc, big.mark = ",", digits = 3), "%")) %>% 
  mutate(val = factor(val, levels = 0:1, labels = c("Unexposed", "Eczema"))) %>% 
  mutate(sleep_class = 
           factor(sleep_class, 
                  levels = c("read", "drugs_poss", "drugs_def", "sleep_def"), 
                  labels = c("All sleep Read codes", "All sleep drugs", "Definite sleep drugs", "Definite sleep drugs + Read codes"))) %>% 
  ggplot(aes(x = val, y = sleep_class, size = Yes, fill = Yes_pc)) +
  geom_point(alpha = 0.6, shape = 21, color = 1) +
  geom_text(aes(label = lab), nudge_y = 0, size = 4, col = 1) +
  scale_size(range = c(5, 20), name = "*n* with code") +
  scale_fill_viridis_c(option = "inferno", name = "*%* with code") +
  labs(x = "Eczema severity", y = "") +
  guides(size = "none") +
  scale_y_discrete(limits=rev) +
  theme_ali() +
  theme(legend.position = "top", 
        legend.title = element_markdown(),
        axis.line = element_line(size = rel(0.5)),
        panel.grid.major.y = element_blank(),
        panel.grid = element_blank())
p3
dev.copy(pdf, here::here("out/supplementary/sleep_data_eczema_yesno.pdf"), width = 8, height = 6); dev.off()


