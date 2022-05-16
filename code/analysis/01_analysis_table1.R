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


if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}
XX <- c("psoriasis", "eczema")
#exposure <- "eczema"

for(exposure in XX) {

ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure)
if (exposure == "eczema") {
  df_anx_split <-
    readRDS(paste0(datapath, "out/df_modelecz_anxiety.rds"))
  df_dep_split <-
    readRDS(paste0(datapath, "out/df_modelecz_depression.rds"))
} else if (exposure == "psoriasis") {
  df_anx_split <-
    readRDS(paste0(datapath, "out/df_modelpso_anxiety.rds"))
  df_dep_split <-
    readRDS(paste0(datapath, "out/df_modelpso_depression.rds"))
}

# psoriasis
# levels(df_anx$severity) <- c("severe", "mild")
# levels(df_dep$severity) <- c("severe", "mild")
table(df_anx_split$severity)
table(df_dep_split$severity)

####
# There seem to be some people with severe records while unexposed controls (before becoming cases). Which seems mental? What is a happening?????
# Currently just suppressing the severe in unexposed but need to investigate
# There also seems to be an absolute mess with the factor levels (from build1.R)

# Build flat baseline char. -----------------------------------------------------------------
build_baseline <- function(df_in = slice(df_anx_split,1:10000)){
  df_in$severity[df_in$exposed == "Unexposed"] <- "None" ## MANUAL SUPPRESSION: not ideal but 94 unexposed~severe in psoriasis (at baseline) so not too big an issue
	if(exposure=="eczema"){
		df_in <- df_in %>% 
			mutate(comorbid = asthma)
						 #severity = factor(severity, levels = c(0, "mild", "severe"), labels=c("Mild", "Moderate", "Severe")))
	}else if(exposure=="psoriasis"){
		df_in <- df_in %>% 
			mutate(comorbid = arthritis)
						 #severity = factor(severity, levels = c(0, "mild"), labels=c("Mild", "Severe")))
	}
	df_out <- df_in %>% 
		arrange(setid, patid, tstart) %>% 
		group_by(setid, patid) %>% 
		slice(1)
	
	tab1 <- df_out %>% 
		ungroup() %>%
		mutate(fup = (enddate-indexdate)/365.25) %>%
		select(setid, patid, exposed, gender, age, agegroup, country, ruc, fup, cal_period, eth_edited, carstairs, ruc, bmi, bmi_cat, obese_cat, alc, smokstatus, smokstatus_nomiss, sleep, sleep_all, comorbid, cci, severity, out) ## need to add DEATH here once it is working properly but seems to have been corrupted by the stsplit (death=1 being copied over multiple lines which is non-sensical)
	table(tab1$exposed, tab1$severity, useNA = "always")
	
	# tab1 <- df_out %>% 
	#   mutate_at("exposed", ~ifelse(. == 0, "Matched controls", paste0("With ", exposure)))
	## reorder value labels
	# tab1 <- tab1 %>% 
	# 	mutate(alc = factor(alc, levels = 0:1, labels = c("No", "Yes")),
	# 				 #smokstatus = factor(smokstatus, levels = 0:1, labels = c("No", "Yes")),
	# 				 sleep = factor(sleep, levels = 0:1, labels = c("No", "Yes")),
	# 				 sleep_all = factor(sleep_all, levels = 0:1, labels = c("No", "Yes")),
	# 				 comorbid = factor(comorbid, levels = 0:1, labels = c("No", "Yes")),
	# 				 cci = factor(cci, levels = 0:2, labels = c("0 Low (0)", "1 Moderate (1-2)", "2 Severe (3 or more)")),
	# 				 out = factor(out, levels = 0:1, labels = c("No", "Yes"))
	# 	)  
	## investigating the unexposed people with severe psoriasis Rx codes
	
	## make table with nice gtsummary 
	table1 <- tab1 %>% 
	  ungroup() %>% 
		select(-patid, -setid, -out) %>%
		tbl_summary(by = exposed,
								statistic = list(all_continuous() ~ "{p50} ({p25}-{p75})",
																 all_categorical() ~ "{n} ({p}%)"),
								digits = all_continuous() ~ 1,
								label = list(gender = "Sex",
														 fup = "Follow-up time (years)",
														 cal_period = "Calendar period",
														 age = "Age",
														 agegroup = "Age (categorised)",
														 eth_edited = "Ethnicity", 
														 bmi = "BMI",
														 bmi_cat = "BMI (categorised)",
														 alc = "Harmful alcohol use", 
														 comorbid = ifelse(exposure=="eczema","Asthma diagnosis","Arthritis diagnosis"),
														 sleep = "Sleep problems", 
														 sleep_all = "Sleep problems (incl. benzos)", 
														 cci = "Charlson's comorbidity index",
														 smokstatus = "Smoking status", 
														 smokstatus_nomiss = "Smoking status (imputed)",
														 carstairs = "Carstairs deprivation quintile",
														 ruc = "Rural/Urban", 
														 country = "Country",
														 severity = paste0("Severe ", exposure) 
								)
		) %>%
		add_overall() %>%
		bold_labels() %>%
		#modify_table_styling(align = "right", columns=6:8) %>%
		modify_footnote(
			all_stat_cols() ~ "Median (IQR) or Frequency (%)"
		) 
	
	table1
}
table1_anx <- build_baseline(df_anx_split)
table1_dep <- build_baseline(df_dep_split)

table1 <- tbl_merge(
	tbls = list(table1_anx, table1_dep),
	tab_spanner = c("**Anxiety cohort**", "**Depression cohort**")
)
table1
table1 %>%
	as_gt() %>%
	gt::gtsave(filename = paste0("tab1_",ABBRVexp,".html"), path = here::here("out/tables"))
}

# Build table 2 - any exposure during follow up ---------------------------
build_followup <- function(df_in){
	#df_in[df_in$exposed == 0, "severity"] <- NA ## MANUAL SUPPRESSION: not ideal but 94 unexposed~severe in psoriasis (at baseline) so not too big an issue
	if(exposure=="eczema"){
		df_in <- df_in %>% 
			mutate(comorbid = asthma)
						 #severity = factor(severity, levels = c(0, "mild", "severe"), labels=c("Mild", "Moderate", "Severe")))
	}else if(exposure=="psoriasis"){
		df_in <- df_in %>% 
			mutate(steroids = NA, # just need a placeholder so the code doesn't break below
						 comorbid = arthritis)
						 #severity = factor(severity, levels = c(0, "mild"), labels=c("Mild", "Severe")))
	}
	df_out <- df_in %>% 
		group_by(setid, patid) %>% 
		mutate(everSev = ifelse(any(severity=="Severe"), 3, 
														ifelse(any(severity == "Moderate"), 2,
																			 ifelse(any(severity == "Mild"), 1, 
																			 					 NA))),
					 everOut = ifelse(any(out==1), 1, 0),
					 everAlc = ifelse(any(alc==1), 1, 0),
					 everSmok = ifelse(any(smokstatus_nomiss == "current smoker"), 3, 
					 											ifelse(any(smokstatus_nomiss == "current or ex-smoker"), 2, 
					 														 ifelse(any(smokstatus_nomiss == "ex-smoker"), 1 ,
					 														 			 ifelse(any(smokstatus_nomiss == "non-smoker"), 0 ,
					 														 			 			 NA)))),
					 everSleep = ifelse(any(sleep==1), 1, 0),
					 everSteroid = ifelse(any(gc90days=="active"), 1, 0),
					 everComorbid = ifelse(any(comorbid==1), 1, 0)
					 ) %>%
		slice(1) %>%
		ungroup() %>%
		mutate(fup = (enddate-indexdate)/365.25) %>%
		select(setid, patid, exposed, everAlc, everSmok, everSleep, everComorbid, everSev, everSteroid, everOut)

	tab2 <- df_out %>% 
		mutate_at("exposed", ~ifelse(. == 0, "Matched controls", paste0("With ", exposure)))
	
	## reorder value labels
	tab2 <- tab2 %>% 
		mutate(everSmok = factor(everSmok, levels = 0:3, labels = c("Non-Smoker", "Ex-Smoker", "Current or ex-smoker", "Current smoker")),
					 everAlc = factor(everAlc, levels = 0:1, labels = c("No", "Yes")),
					 everSleep = factor(everSleep, levels = 0:1, labels = c("No", "Yes")),
					 everComorbid = factor(everComorbid, levels = 0:1, labels = c("No", "Yes")),
					 everSev = factor(everSev, levels = 1:3, labels = c("Mild", "Moderate","Severe")),
					 everSteroid = factor(everSteroid, levels = 0:1, labels = c("None", "At least 1")),
					 everOut = factor(everOut, levels = 0:1, labels = c("No", "Yes"))
					 )
	
	## investigating the unexposed people with severe psoriasis Rx codes
	table(tab2$exposed, tab2$everSev)
	
	## make table with nice gtsummary 
	if(exposure=="psoriasis"){tab2 <- select(tab2, -everSteroid)}
	table2 <- tab2 %>% 
		select(-patid, -setid, -everOut) %>%
		tbl_summary(by = exposed,
								statistic = list(all_continuous() ~ "{p50} ({p25}-{p75})",
																 all_categorical() ~ "{n} ({p}%)"),
								digits = all_continuous() ~ 1,
								label = list(gender = "Sex",
														 fup = "Follow-up time (years)",
														 cci = "Charlson's comorbidity index",
														 age = "Age",
														 agegroup = "Age (categorised)",
														 eth_edited = "Ethnicity", 
														 bmi = "BMI",
														 bmi_cat = "BMI (categorised)",
														 #everOut = outcome_cohort,
														 everAlc = "Harmful alcohol use", 
														 everComorbid = ifelse(exposure=="eczema","Asthma diagnosis","Arthritis diagnosis"),
														 everSleep = "Sleep problems", 
														 everSmok = "Smoking status", 
														 carstairs = "Carstairs deprivation quintile",
														 ruc = "Rural/Urban", 
														 everSteroid = "High dose oral steroid prescription",
														 everSev = paste0("Severe ", exposure)
														 )
		) %>%
		add_overall() %>%
		bold_labels() %>%
		modify_table_styling(align = "right", columns=6:8) %>%
		modify_footnote(
			all_stat_cols() ~ "Median (IQR) or Frequency (%)"
		) 
	
	table2
}
table2_anx <- build_followup(df_anx_split)
table2_dep <- build_followup(df_dep_split)

table2 <- tbl_merge(
		tbls = list(table2_anx, table2_dep),
		tab_spanner = c("**Anxiety cohort**", "**Depression cohort**")
	)
table2
table2 %>%
	as_flex_table() %>%
	flextable::save_as_docx(path = here::here("out/tables", paste0("tab2_",ABBRVexp,"_fup.docx")))
table2 %>%
	as_gt() %>%
	gt::gtsave(filename = paste0("tab2_",ABBRVexp,"_fup.html"), path = here::here("out/tables"))


# Build table 3 - person years ---------------------------------------------------
## will need this little function later to summ follow up by factor vars
summ_pyars <- function(V = "alc", df = tab1){
	df %>%
		group_by(exposed, get(V)) %>%
		summarise(pyar = sum(time)) %>%
		rename(var = `get(V)`) %>%
		mutate(name = V,
					 group_pyar = sum(pyar)) %>%
		ungroup()
}

## and this is the main function 
build_pyears <- function(df_in){
	#df_in[df_in$exposed == 0, "severity"] <- NA ## MANUAL SUPPRESSION: not ideal but 94 unexposed~severe in psoriasis (at baseline) so not too big an issue
	if(exposure=="eczema"){
		df_in <- df_in %>% 
			mutate(comorbid = asthma)
						 #severity = factor(severity, levels = c(0, "mild", "severe"), labels=c("Mild", "Moderate", "Severe")))
	}else if(exposure=="psoriasis"){
		df_in <- df_in %>% 
			mutate(steroids = NA, # just need a placeholder so the code doesn't break below
						 comorbid = arthritis)
						 #severity = factor(severity, levels = c(0, "mild"), labels=c("Mild", "Severe")))
	}

	tab1 <- df_in %>% 
		as_tibble() %>% 
		ungroup() %>%
		mutate_at("exposed", ~ifelse(. == 0, "Matched controls", "With exposure")) %>%
		arrange(setid, patid) %>%
		mutate(time = as.numeric(tstop-tstart)/365.25) %>%
		mutate(sleep = as.factor(sleep),
					 country = as.factor(country),
					 gc90day = as.factor(gc90days))

	#summ_pyars("comorbid")
	
	x <- tab1 %>%
						select_if(is.factor) %>%
						select(-arthritis, -asthma, -age_cat) %>%
						names() %>% as.list()
	
	pyars_table <- map(x, summ_pyars, df = tab1) %>% ## to use the little function defined above
		do.call(rbind, .)
	
	gt_pyars_table <- pyars_table %>%
		group_by(exposed) %>%
		mutate(pc_pyar = (pyar/group_pyar)*100) %>%
		pivot_wider(
			id_cols = c(name, var), 
			names_from = c(exposed), 
			values_from = c(pyar, pc_pyar)) %>%
		clean_names() %>%
		mutate_at("var", ~stringr::str_to_title(.)) %>%
		mutate(
			new_lab = case_when(
			name == "gender" ~ "Sex",
			name == "agegroup" ~ "Age (categorised)",
			name == "eth_edited" ~ "Ethnicity", 
			name == "bmi_cat" ~ "BMI (categorised)",
			name == "alc" ~ "Harmful alcohol use", 
			name == "comorbid" ~ ifelse(exposure=="eczema","Asthma diagnosis","Arthritis diagnosis"),
			name == "cci" ~ "Charlson's comorbidity index",
			name == "sleep" ~ "Sleep problems", 
			name == "smokstatus_nomiss" ~ "Smoking status", 
			name == "carstairs" ~ "Carstairs deprivation quintile",
			name == "cal_period" ~ "Calendar period", 
			name == "gc90day" ~ "Oral steroid prescription", 
			name == "ruc" ~ "Rural/Urban", 
			name == "country" ~ "Country", 
			name == "severity" ~ "Severity", 
			)
			) %>%
		select(-name) %>%
		arrange(new_lab)
	
	gt_pyars_table %>%
		gt(
			rowname_col = "var",
			groupname_col = "new_lab"
		) %>%
		tab_stubhead(
			label = "Characteristic"
		) %>%
		fmt_number(
			columns = where(is.numeric),
			decimals = 1, 
			use_seps = T
		) %>%
		cols_merge(
			columns = contains("matched_controls"),
			pattern = "{1} ({2}%)"
		) %>%
		cols_merge(
			columns = contains("with_exposure"),
			pattern = "{1} ({2}%)"
		) %>%
		tab_style(
			style = cell_text(weight = "bold"),
			locations = cells_row_groups(groups = everything())
		) %>%
		cols_label(
			pyar_matched_controls = md("**Matched controls**"),
			pyar_with_exposure = md(paste0("**With ",exposure,"**"))
		)
}
table3_anx <- build_pyears(df_in=df_anx_split)
table3_dep <- build_pyears(df_in=df_dep_split)

table3_anx %>%
	gt::gtsave(filename = paste0("tab3_",ABBRVexp,"_pyars_anxiety.html"), path = here::here("out/tables")) 

table3_dep %>%
	gt::gtsave(filename = paste0("tab3_",ABBRVexp,"_pyars_depression.html"), path = here::here("out/tables")) 

}