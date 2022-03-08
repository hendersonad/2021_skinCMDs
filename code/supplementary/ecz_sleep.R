# install.packages("rtools")
# install.packages("gt")
# install.packages("gtsummary")
# install.packages("janitor")
# install.packages("here")

library(tidyverse)
library(here)
library(flextable)
library(gt)
library(gtsummary)
library(janitor)
library(timetk)
library(skimr)
library(glue)
library(gridExtra)


if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		#datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

exposure <- "eczema"
ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure)

df_anx_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-anxiety_split.rds"))
df_dep_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_split.rds"))
cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))

# merge on sleep codes to get Source of sleep diagnosis -------------------

ecz_sleep_def <- haven::read_dta(paste0(datapath, "out/variables-ecz-sleep-definite.dta"))
ecz_sleep_all <- haven::read_dta(paste0(datapath, "out/variables-ecz-sleep-all.dta"))

duplicated(ecz_sleep_def$patid) %>% sum() # all unique so that's good

add_text <- function(x, text){
  paste0(x,"_",text)
}
merge_ecz_sleep_def <- ecz_sleep_def %>% 
  select(patid, ingredient, src, readterm) %>% 
  rename_with(add_text, text = "def")
merge_ecz_sleep_all <- ecz_sleep_all %>% 
  select(patid, ingredient, src, readterm) %>% 
  rename_with(add_text, text = "all")

df_anx_split_withsrc <- df_anx_split %>% 
  ungroup() %>% 
  select(setid, patid, tstart, exposed, severity, sleep, sleep_all) %>% 
  left_join(merge_ecz_sleep_def, by = c("patid" = "patid_def")) %>% 
  left_join(merge_ecz_sleep_all, by = c("patid" = "patid_all"))  
  
df_dep_split_withsrc <- df_dep_split %>% 
  ungroup() %>% 
  select(setid, patid, tstart, exposed, severity, sleep, sleep_all) %>% 
  left_join(merge_ecz_sleep_def, by = c("patid" = "patid_def")) %>% 
  left_join(merge_ecz_sleep_all, by = c("patid" = "patid_all"))  

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

for (outcome in c("depression", "anxiety")) {
	if(outcome == "depression"){
		df_in = df_dep_split_withsrc
	}
	if(outcome == "anxiety"){
		df_in = df_anx_split_withsrc
	}
	df_out <- df_in %>% 
		arrange(setid, patid, tstart) %>% 
		group_by(setid, patid) %>% 
		slice(1)
	
	df_out$src_def[is.na(df_out$src_def)] <- 0
	df_out$src_all[is.na(df_out$src_all)] <- 0
	
  # all sleep codes - drugs only --------------------------------------------
	df_out_drugs_all <- df_out %>% 
	  ungroup() %>% 
	  filter(src_all != 1)
	tab_sleep_drugs_all <- twoXtwo(df = df_out_drugs_all, exp = "severity", out = "sleep_all") %>% drop_na() %>% mutate(sleep_class = "drugs")
	
  # all sleep codes - read only --------------------------------------------
	df_out_read_all <- df_out %>% 
	  ungroup() %>% 
	  filter(src_all != 2)
	tab_sleep_read_all <- twoXtwo(df = df_out_read_all, exp = "severity", out = "sleep_all") %>% drop_na() %>% mutate(sleep_class = "read code")
	
	# def sleep codes - drugs only --------------------------------------------
	df_out_drugs_def <- df_out %>% 
	  ungroup() %>% 
	  filter(src_def != 1)
	tab_sleep_drugs_def <- twoXtwo(df = df_out_drugs_def, exp = "severity", out = "sleep") %>% drop_na() %>% mutate(sleep_class = "drugs")
	
  # def sleep codes - read only --------------------------------------------
	df_out_read_def <- df_out %>% 
	  ungroup() %>% 
	  filter(src_def != 2)
	tab_sleep_read_def <- twoXtwo(df = df_out_read_def, exp = "severity", out = "sleep") %>% drop_na() %>% mutate(sleep_class = "read code")
	
	tab_sleep_out <- bind_rows(tab_sleep_drugs_def, tab_sleep_read_def,
	                               tab_sleep_drugs_all, tab_sleep_read_all)
	
	gt_sleep <- tab_sleep_out %>%
		select(-exposure,-Miss, -sleep_class) %>%
		gt::gt() %>%
		tab_row_group(label = "Definite sleep drugs",
									rows = 1:4) %>%
		tab_row_group(label = "Definite sleep Read codes",
									rows = 5:8) %>%
		tab_row_group(label = "All sleep drugs (incl. benzo)",
									rows = 9:12) %>%
		tab_row_group(label = "All sleep Read codes",
									rows = 13:16) %>%
		gt::tab_header(title = "Sleep problems by eczema severity") %>%
		gt::fmt_number(columns = c(3, 5), decimals = 1) %>%
		gt::fmt_number(columns = c(2, 4), decimals = 0) %>%
		gt::data_color(
			columns = c(Yes_pc),
			colors = scales::col_numeric(
				palette = paletteer::paletteer_c(palette = "harrypotter::gryffindor",
																				 n = 100) %>% as.character(),
				domain = c(
					min(tab_sleep_out$Yes_pc, na.rm = T),
					max(tab_sleep_out$Yes_pc, na.rm = T)
				)
			)
		) %>%
		cols_label(val = "Severity",
							 No_pc = "%",
							 Yes_pc = "%") 
	
	
	
	gt_sleep %>% 
		gt::gtsave(
			filename =  paste0("eczema_sleep_", outcome,"_v2.html"),
			path = here::here("out/supplementary")
		)
}

# what ingredients and diagnoses are most prevalent  ----------------------

cohort_join <- cohort %>% 
  select(setid, patid, exposed)
exposed_tab <- table(cohort_join$exposed)
n_exp <- exposed_tab[2]
n_un <- exposed_tab[1]

ecz_sleep_def_describe <- merge_ecz_sleep_def %>% 
  mutate(desc = ifelse(ingredient_def == "", readterm_def, ingredient_def)) %>% 
  mutate(desc = str_remove(desc, "\\[d\\]")) %>% 
  mutate(desc = str_remove(desc, "c/o - ")) %>% 
  mutate(desc = str_remove(desc, "\\[x\\]"))  
  
ecz_sleep_def_describe2 <- ecz_sleep_def_describe %>% 
  left_join(cohort_join, by = c("patid_def" = "patid"))

ecz_sleep_def_plot <- ecz_sleep_def_describe2 %>%
  group_by(desc, src_def, exposed) %>% 
  summarise(n = n()) %>% 
  filter(n>1000) %>% 
  arrange(-n) %>% 
  mutate(rank = fct_reorder(as.factor(desc), n)) 

ecz_sleep_def_plot$rank <- fct_reorder(ecz_sleep_def_plot$rank, ecz_sleep_def_plot$n)

ecz_sleep_def_plot$src <- as.factor(ecz_sleep_def_plot$src_def)
levels(ecz_sleep_def_plot$src) <- c("Read code", "Prescription")
ecz_sleep_def_plot$exposed <- as.factor(ecz_sleep_def_plot$exposed)
levels(ecz_sleep_def_plot$exposed) <- c("Unexposed", "With eczema")

ecz_sleep_def_plot$bigN <- ifelse(ecz_sleep_def_plot$exposed=="Unexposed", n_un, n_exp)
ecz_sleep_def_plot$pc <- round((ecz_sleep_def_plot$n/ecz_sleep_def_plot$bigN)*100, 1)
ecz_sleep_def_plot$prettyN <- paste0(prettyNum(ecz_sleep_def_plot$n, big.mark = ","), " (", ecz_sleep_def_plot$pc, "%)")
ecz_sleep_def_plot %>% 
  group_by(exposed) %>% 
  summarise(sum(pc))
ggplot(ecz_sleep_def_plot, aes(x = rank, y = n, group = src, col = src, fill = src)) +
  geom_col(alpha = 0.4) + 
  geom_text(aes(y = 2.1e5, label = prettyN), hjust = "right") +
  coord_flip() +
  facet_wrap(~exposed, ncol = 2) + 
  labs(x = "Sleep code") + 
  scale_y_continuous(name="Number of instances", labels = scales::comma) + 
  theme_ali() +
  theme(legend.position = "bottom",
        strip.background = element_blank())

dev.copy(pdf, here::here("out/supplementary/sleep_definite.pdf"), width = 8, height = 6); dev.off()

ecz_sleep_all_describe <- merge_ecz_sleep_all %>% 
  mutate(desc = ifelse(ingredient_all == "", readterm_all, ingredient_all)) %>% 
  mutate(desc = str_remove(desc, "\\[d\\]")) %>% 
  mutate(desc = str_remove(desc, "c/o - ")) %>% 
  mutate(desc = str_remove(desc, "\\[x\\]"))  

ecz_sleep_all_describe2 <- ecz_sleep_all_describe %>% 
  left_join(cohort_join, by = c("patid_all" = "patid"))

ecz_sleep_all_plot <- ecz_sleep_all_describe2 %>%
  group_by(desc, src_all, exposed) %>% 
  summarise(n = n()) %>% 
  filter(n>1000) %>% 
  arrange(-n) %>% 
  mutate(rank = fct_reorder(as.factor(desc), n)) 

ecz_sleep_all_plot$rank <- fct_reorder(ecz_sleep_all_plot$rank, ecz_sleep_all_plot$n)

ecz_sleep_all_plot$src <- as.factor(ecz_sleep_all_plot$src_all)
levels(ecz_sleep_all_plot$src) <- c("Read code", "Prescription")
ecz_sleep_all_plot$exposed <- as.factor(ecz_sleep_all_plot$exposed)
levels(ecz_sleep_all_plot$exposed) <- c("Unexposed", "With eczema")

ecz_sleep_all_plot$prettyN <- prettyNum(ecz_sleep_all_plot$n, big.mark = ",")

ecz_sleep_all_plot$bigN <- ifelse(ecz_sleep_all_plot$exposed=="Unexposed", n_un, n_exp)
ecz_sleep_all_plot$pc <- round((ecz_sleep_all_plot$n/ecz_sleep_all_plot$bigN)*100, 1)
ecz_sleep_all_plot$prettyN <- paste0(prettyNum(ecz_sleep_all_plot$n, big.mark = ","), " (", ecz_sleep_all_plot$pc, "%)")

ggplot(ecz_sleep_all_plot, aes(x = rank, y = n, group = src, col = src, fill = src)) +
  geom_col(alpha = 0.4) + 
  geom_text(aes(y = 4e5, label = prettyN), hjust = "right") +
  coord_flip() +
  facet_wrap(~exposed, ncol = 2) + 
  labs(x = "Sleep code") + 
  scale_y_continuous(name="Number of instances", labels = scales::comma) + 
  theme_ali() +
  theme(legend.position = "bottom",
        strip.background = element_blank())

dev.copy(pdf, here::here("out/supplementary/sleep_all.pdf"), width = 8, height = 6); dev.off()

