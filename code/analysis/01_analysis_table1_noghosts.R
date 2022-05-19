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
  df_anx_split <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_anxiety.rds"))
  df_dep_split <- readRDS(paste0(datapath, "out/df_model", ABBRVexp, "_depression.rds"))
  var_consbeforeindex <-
    haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-consultations-yrbeforeindex-3yrs.dta"))
  
  ####
  # There seem to be some people with severe records while unexposed controls (before becoming cases). Which seems mental? What is a happening?????
  # Currently just suppressing the severe in unexposed but need to investigate
  # There also seems to be an absolute mess with the factor levels (from build1.R)

  # Build flat baseline char. -----------------------------------------------------------------
  build_baseline <- function(df_in = slice(df_anx_split,1:1000)){
  	df_in[df_in$exposed == 0, "severity"] <- "none" ## MANUAL SUPPRESSION: not ideal but 94 unexposed~severe in psoriasis (at baseline) so not too big an issue
  	if(exposure=="eczema"){
  		df_in <- df_in %>% 
  			mutate(comorbid = asthma)
  	}else if(exposure=="psoriasis"){
  		df_in <- df_in %>% 
  			mutate(comorbid = arthritis)
    }
  	
  	# add in filtering for consultation <= 1year before index -----------------
  	df_model_noghosts <- df_in %>% 
  	  left_join(dplyr::select(var_consbeforeindex, -indexdate), by = "patid") %>% 
  	  mutate(pre_cons = replace_na(consyrbeforeindex, 0)) %>% 
  	  filter(pre_cons == 1)
  	
  	# only keep valid sets ----------------------------------------------------
  	validsets_noghosts <- df_model_noghosts %>% 
  	  group_by(setid) %>% 
  	  summarise(mean = mean(exposed == str_to_title(exposure))) %>% 
  	  filter(mean > 0 & mean < 1)
  	validsets <- validsets_noghosts$setid
  	
  	df_validsets_noghosts <- df_model_noghosts %>% 
  	  filter(setid %in% validsets)
  	
  	df_out <- df_validsets_noghosts %>% 
  		arrange(setid, patid, tstart) %>% 
  		group_by(setid, patid) %>% 
  		slice(1)
  	
  	tab1 <- df_out %>%
  	  ungroup() %>%
  	  mutate(fup = (enddate - indexdate) / 365.25) %>%
  	  select(setid,patid,exposed,gender,age,agegroup,country,ruc,fup,cal_period,eth_edited,carstairs,ruc,bmi,bmi2,bmi_cat,alc,smokstatus,sleep,sleep_all,comorbid,cci,severity,out) ## need to add DEATH here once it is working properly but seems to have been corrupted by the stsplit (death=1 being copied over multiple lines which is non-sensical)
  	table(tab1$exposed, tab1$severity, useNA = "always")
  	
  	## make table with nice gtsummary 
  	table1 <- tab1 %>%
  	  ungroup() %>%
  	  select(-patid,-setid,-out) %>%
  	  tbl_summary(
  	    by = exposed,
  	    statistic = list(
  	      all_continuous() ~ "{p50} ({p25}-{p75})",
  	      all_categorical() ~ "{n} ({p}%)"
  	    ),
  	    digits = all_continuous() ~ 1,
  	    label = list(
  	      gender = "Sex",
  	      fup = "Follow-up time (years)",
  	      cal_period = "Calendar period",
  	      age = "Age",
  	      agegroup = "Age (categorised)",
  	      eth_edited = "Ethnicity",
  	      bmi = "BMI",
  	      bmi2 = "BMI (centred)",
  	      bmi_cat = "BMI (categorised)",
  	      alc = "Harmful alcohol use",
  	      comorbid = ifelse(
  	        exposure == "eczema",
  	        "Asthma diagnosis",
  	        "Arthritis diagnosis"
  	      ),
  	      sleep = "Sleep problems",
  	      sleep_all = "Sleep problems (incl. benzos)",
  	      cci = "Charlson's comorbidity index",
  	      smokstatus = "Smoking status",
  	      carstairs = "Carstairs deprivation quintile",
  	      ruc = "Rural/Urban",
  	      country = "Country",
  	      severity = paste0("Severe ", exposure)
  	    )
  	  ) %>%
  	  add_overall() %>%
  	  bold_labels() %>%
  	  modify_footnote(all_stat_cols() ~ "Median (IQR) or Frequency (%)")
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
  	gt::gtsave(filename = paste0("tab1_",ABBRVexp,"_noghosts.html"), path = here::here("out/tables"))
}
