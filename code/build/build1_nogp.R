library(here, haven, tidyverse)
library(survival, survminer)
library(gt, gtsummary)
library(janitor, flextable)
library(skimr)
library(clipr)
library(KMunicate)
library(data.table)
library(arrow)
library(formattable)

if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

XX <- c("psoriasis")

for(exposure in XX){
	#exposure <- XX[1]
	ABBRVexp <- str_sub(exposure,1 ,3)
	.dib(exposure) 

	# import data -------------------------------------------------------------
	patient <- haven::read_dta(paste0(datapath, "in/nogp/Patient_extract_", ABBRVexp, "_extract3_nogp_1.dta"))
	cohort <- haven::read_dta(paste0(datapath, "out/nogp/getmatchedcohort-", exposure, "-main-mhealth-nogp.dta"))
	var_depressionDef <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-depression-definite-nogp.dta"))
	var_anxietyDef <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-anxiety-definite-nogp.dta"))
	var_depressionAll <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-depression-all-nogp.dta"))
	var_anxietyAll <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-anxiety-all-nogp.dta"))
	var_depressionCensor <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-depression-censor-nogp.dta"))
	var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-anxiety-censor-nogp.dta"))
	var_agesexgp <- haven::read_dta(paste0(datapath, "out/nogp/variables-", ABBRVexp, "-age-sex-nogp", ".dta"))
			
			# censor outcome ----------------------------------------------------------
			var_dep <- var_depressionAll %>%
			  rename(depdate = eventdate) %>%
			  left_join(var_depressionCensor, by = "patid") %>%
			  rename(censorDepDate = eventdate) %>%
			  mutate_at("censorDepDate", ~ifelse(is.na(.), as.Date("3000-01-01"), as.Date(.))) %>%
			  mutate_at("censorDepDate", ~as.Date(., origin = "1970-01-01")) %>%
			  mutate(dep = case_when(
			    depdate < censorDepDate &
			      !is.na(depdate) ~ 1,
			    TRUE ~ 0
			  )) %>%
			  select(patid, case = dep, eventdate = depdate, censorDepDate, readterm=readterm.x)
			
			depression_censoring <- var_dep %>% 
			  filter(censorDepDate == eventdate)
			
			var_anx <- var_anxietyAll %>%
			  rename(anxdate = eventdate) %>%
			  left_join(var_anxietyCensor, by = "patid") %>%
			  rename(censorAnxDate = eventdate) %>%
			  mutate_at("censorAnxDate", ~ifelse(is.na(.), as.Date("3000-01-01"), as.Date(.))) %>%
			  mutate_at("censorAnxDate", ~as.Date(., origin = "1970-01-01")) %>%
			  mutate(anx = case_when(
			    anxdate < censorAnxDate &
			      !is.na(anxdate) ~ 1,
			    TRUE ~ 0
			  )) %>%
			  select(patid, case = anx, eventdate = anxdate, censorAnxDate, readterm=readterm.x)
			
			#rm(var_depressionAll, var_depressionCensor,var_anxietyAll, var_anxietyCensor)
			
			ggplot(var_anx, aes(x=eventdate, group = case)) +
			  geom_histogram(bins = 50, fill="gray80", col="gray30") +
			  facet_wrap(~case)+
			  theme_bw()
			ggplot(var_dep, aes(x=eventdate, group = case)) +
			  geom_histogram(bins = 50, fill="gray80", col="gray30") +
			  facet_wrap(~case)+
			  theme_bw()
		  
			cohort_edited <- cohort
		  cohort_edited <- cohort_edited %>%
		    select(setid, patid, exposed, indexdate, enddate) %>%
		    left_join(var_agesexgp, by = "patid")
			
			# ANXIETY -----------------------------------------------------------------
			## make into tmerge format
			select_vars <- c("enddate", "eventdate")
			
			### First got to censor all the sets 
			an_anxiety <- cohort_edited %>% 
			  ungroup() %>%
			  left_join(var_anx, by = "patid") %>%
			  mutate(newenddate = apply(select(., all_of(select_vars)), 1, FUN = min, na.rm = T)) %>%	
			  filter(newenddate>indexdate) %>%
			  select(-newenddate) ## have used this data enough now
			
			## check valid sets
			.dib(paste0("no of unique sets BEFORE filtering: ", length(unique(cohort_edited$setid))))
			invalid_sets <- an_anxiety %>%
			  group_by(setid) %>%
			  summarise(valid=max(exposed)) %>%
			  filter(valid==0) %>%
			  pull(setid)
			an_anxiety <- an_anxiety %>%
			  filter(!setid %in% invalid_sets)
			.dib(paste0("no of unique sets AFTER filtering: ", length(unique(an_anxiety$setid))))
			.dib(paste0("no of sets REMOVED: ", length(unique(cohort_edited$setid))-length(unique(an_anxiety$setid))))
			
			## Store the edited anxiety dataset 
			# and tidy (get rid of anxiety data to merge back on in a minute using tmerge)
			an_anxiety <- an_anxiety %>% 
			  select(setid, patid, exposed, indexdate, enddate, realyob, gender, pracid)
			if(export_datasets) {
			  saveRDS(an_anxiety,
			          file = paste0(datapath, "out/nogp/", ABBRVexp, "_an_anxiety.rds"))
			}
			if(sum(grepl(pattern = "cohort_edited", x = ls()))==0){
			  an_anxiety <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "_an_anxiety.rds"))
			}
			
			# DEPRESSION --------------------------------------------------------------
			## make into tmerge format
			select_vars <- c("enddate", "eventdate")
			
			an_depression <- cohort_edited %>% 
			  ungroup() %>%
			  left_join(var_dep, by = "patid") %>%
			  mutate(newenddate = apply(select(., all_of(select_vars)), 1, FUN = min, na.rm = T)) %>%
			  #mutate(newenddate = pmin(enddate, depdate, na.rm=T)) %>%
			  filter(newenddate>indexdate) %>%
			  select(-newenddate) ## have used this data enough now
			
			## check valid sets
			.dib(paste0("no of unique sets BEFORE filtering: ", length(unique(cohort_edited$setid))))
			invalid_sets <- an_depression %>%
			  group_by(setid) %>%
			  summarise(valid=max(exposed)) %>%
			  filter(valid==0) %>%
			  pull(setid)
			an_depression <- an_depression %>%
			  filter(!setid %in% invalid_sets)
			.dib(paste0("no of unique sets AFTER filtering: ", length(unique(an_depression$setid))))
			.dib(paste0("no of sets REMOVED: ", length(unique(cohort_edited$setid))-length(unique(an_depression$setid))))
			
			## Store the edited depression dataset 
			# and tidy (get rid of anxiety data to merge back on in a minute using tmerge)
			an_depression <- an_depression %>% 
			  select(setid, patid, exposed, indexdate, enddate, realyob, gender, pracid)
			
			if(export_datasets) {
			  saveRDS(an_depression,
			          file = paste0(datapath, "out/nogp/", ABBRVexp, "_an_depression.rds"))
			}
			
			if(sum(grepl(pattern = "cohort_edited", x = ls()))==0){
			  an_depression <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "_an_depression.rds"))
			}
			
			
			# tmerge function ---------------------------------------------------------
			build_tmerge <- function(exp = 1, df_outcome = an_anxiety, outcome = var_anx){
			  df_exp <- df_outcome %>%
			    filter(exposed==exp)
			  
			  out_tmerge0 <- tmerge(df_exp, df_exp, id=patid, tstart = indexdate, tstop = enddate) #This first step sets the range of follow up
			  out_tmerge1 <- tmerge(out_tmerge0, outcome, id=patid, out = event(eventdate))
			  
			  ## Restrict to first occurence of outcome ## post-processing to keep only records upto event if have the event. 
			  out_tmerge9 <- out_tmerge1 %>%
			    group_by(patid) %>%
			    slice(if(any(out == 1)){1:which.max(out == 1)}else{row_number()}) %>% 
			    ungroup()
			  out_tmerge9
			}
			anxiety_exposed   <- build_tmerge(exp=1, df_outcome = an_anxiety, outcome = var_anx)
			anxiety_unexposed <- build_tmerge(exp=0, df_outcome = an_anxiety, outcome = var_anx)
			anxiety_full <- anxiety_exposed %>%
			  bind_rows(anxiety_unexposed)
			dim(anxiety_full)
			
			depression_exposed   <- build_tmerge(exp=1, df_outcome = an_depression, outcome = var_dep)
			depression_unexposed <- build_tmerge(exp=0, df_outcome = an_depression, outcome = var_dep)
			depression_full <- depression_exposed %>%
			  bind_rows(depression_unexposed)
			dim(depression_full)
			
			if(export_datasets){
			  saveRDS(anxiety_full, file = paste0(datapath, "out/nogp/", ABBRVexp, "-anxiety_full.rds")) ## some temp stats from wrong run of tmerge commands - 3541262 (16 cols). Now 7234333 (because I forgot to change exp=0 in the unexposed run like a fool)
			  saveRDS(depression_full, file = paste0(datapath, "out/nogp/", ABBRVexp, "-depression_full.rds")) ## some temp stats from wrong run of tmerge commands - 6286583 (16 cols)		
			}	
			
			## age as underlying timescale
			anxiety_full$dob <- as.Date(paste(anxiety_full$realyob, "06-01", sep="-"))
			anxiety_full$age_at_index <- as.numeric(anxiety_full$indexdate - anxiety_full$dob)/365.25
			anxiety_full$age_at_end <- as.numeric(as.Date(anxiety_full$enddate) - anxiety_full$dob)/365.25
			anxiety_full$age_at_tstart <- as.numeric(anxiety_full$tstart - anxiety_full$dob)/365.25
			anxiety_full$age_at_tstop <- as.numeric(as.Date(anxiety_full$tstop) - anxiety_full$dob)/365.25
			anxiety_full$t <- anxiety_full$age_at_tstop - anxiety_full$age_at_tstart
			anxiety_full$anx <- anxiety_full$out
			## age as underlying timescale
			depression_full$dob <- as.Date(paste(depression_full$realyob, "06-01", sep="-"))
			depression_full$age_at_index <- as.numeric(depression_full$indexdate - depression_full$dob)/365.25
			depression_full$age_at_end <- as.numeric(as.Date(depression_full$enddate) - depression_full$dob)/365.25
			depression_full$age_at_tstart <- as.numeric(depression_full$tstart - depression_full$dob)/365.25
			depression_full$age_at_tstop <- as.numeric(as.Date(depression_full$tstop) - depression_full$dob)/365.25
			depression_full$t <- depression_full$age_at_tstop - depression_full$age_at_tstart
			depression_full$dep <- depression_full$out
			
			# merge on static vars ----------------------------------------------------
			rm(an_anxiety, an_depression, anxiety_exposed, anxiety_unexposed, depression_exposed, depression_unexposed,
			   cohort, cohort_edited)
			
			if(sum(grepl(pattern = "anxiety_full", x = ls()))==0){
			  anxiety_full <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "-anxiety_full.rds"))
			  depression_full <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "-depression_full.rds"))
			  anxiety_full$dob <- as.Date(paste(anxiety_full$realyob, "06-01", sep="-"))
			  depression_full$dob <- as.Date(paste(depression_full$realyob, "06-01", sep="-"))
			}
			
			merge_static <- function(df_out = depression_full){
			    df_static <- df_out %>% 
			      select(setid, patid, pracid, exposed, indexdate, enddate, dob, gender, tstart, tstop, out) 
			  
			  df_static <- df_static %>% 
			    mutate(
			      age=as.numeric(tstart-dob)/365.25,
			      age_cat=case_when((as.numeric(tstart)-as.numeric(dob))/365.25 < 30 ~ "18 - 29",
			                        (as.numeric(tstart)-as.numeric(dob))/365.25 < 40 ~ "30 - 39",
			                        (as.numeric(tstart)-as.numeric(dob))/365.25 < 50 ~ "40 - 49",
			                        (as.numeric(tstart)-as.numeric(dob))/365.25 < 66 ~ "50 - 65",
			                        (as.numeric(tstart)-as.numeric(dob))/365.25 >= 66 ~ "66+"),
			      age_cat=factor(age_cat, levels = c("18 - 29", 
			                                         "30 - 39", 
			                                         "40 - 49", 
			                                         "50 - 65", 
			                                         "66+"))
			    )
			  
			  df_static
			}
			anxiety_static <- merge_static(anxiety_full)
			depression_static <- merge_static(depression_full)
			
		  saveRDS(anxiety_static, file = paste0(datapath, "out/nogp/", ABBRVexp, "-anxiety_static.rds")) 
		  saveRDS(depression_static, file = paste0(datapath, "out/nogp/", ABBRVexp, "-depression_static.rds")) 
			
			build_split <- function(df_in = anxiety_static){
			  tosplit <- df_in %>% 
			    arrange(setid, patid)
			  
			  df_split_cal <- survSplit(Surv(time = as.numeric(tstart),
			                                 time2 = as.numeric(tstop),
			                                 event = out) ~ .,
			                            data=tosplit, 
			                            cut=as.numeric(as.Date(c("2002-01-01","2008-01-01","2014-01-01"))),
			                            episode="cal_period")
			  
			  ## store dob before it gets messed up by split
			  df_split_cal <- df_split_cal %>%
			    mutate(dobNum = as.numeric(dob))
			  ## split on age
			  df_split_age <- survSplit(Surv(time = as.numeric(tstart), 
			                                 time2 = as.numeric(tstop),
			                                 origin = dobNum,
			                                 event = out) ~ ., 
			                            data=df_split_cal, 
			                            cut=c(40, 50, 60, 70, 80)*365.25,
			                            episode="agegroup")
			  #Make factors
			  df_split_age$cal_period <- factor(df_split_age$cal_period,levels = 1:4,
			                                    labels = c("1997-2002", "2002-08", "2008-14", "2014-19"))
			  df_split_age$agegroup <- factor(df_split_age$agegroup,levels = c(1, 2, 3, 4, 5, 6),
			                                  labels = c("18-39", "40-49", "50-59", "60-69", "70-79", "80+"))
			  
			  df_split_age
			}
			rm(anxiety_full, depression_full, var_agesexgp)
			
			anxiety_split <- build_split(anxiety_static)
			depression_split <- build_split(depression_static)
			
		  saveRDS(anxiety_split, file = paste0(datapath, "out/nogp/", ABBRVexp, "-anxiety_split.rds")) 
		  saveRDS(depression_split, file = paste0(datapath, "out/nogp/", ABBRVexp, "-depression_split.rds")) 
}


for(exposure in XX){
  #exposure <- XX[1]
  ABBRVexp <- str_sub(exposure,1 ,3)
  .dib(exposure) 
  
  anxiety_split <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "-anxiety_split.rds"))
  depression_split <- readRDS(paste0(datapath, "out/nogp/", ABBRVexp, "-depression_split.rds"))
  for(outcome in c("depression", "anxiety")){
    if (outcome == "anxiety") {
      df_model <- anxiety_split
    } else if (outcome == "depression") {
      df_model <- depression_split
    }
    
    # Bit of variable formatting for output in regression tables --------------
    df_model$t <-
      as.numeric(df_model$tstop - df_model$tstart)
    
    df_model$gc90days <-
      factor(df_model$gc90days, levels = c("not active", "active"))
    
    df_model$bmi2 <-
      df_model$bmi - mean(df_model$bmi, na.rm = T)
    
    df_model <- df_model %>% 
      mutate(exposed = case_when(
        exposed == 0 ~ "Unexposed",
        exposed == 1 ~ stringr::str_to_title(paste0(exposure))
      ))
    df_model$exposed <- factor(df_model$exposed, levels = c("Unexposed", str_to_title(exposure)))
    var_label(df_model$exposed) <- "Exposure"
    var_label(df_model$carstairs) <- "Carstairs index of deprivation"
    levels(df_model$carstairs)[1] <- "1 (least deprived)"
    levels(df_model$carstairs)[5] <- "5 (most deprived)"
    var_label(df_model$cal_period) <- "Calendar Period"
    var_label(df_model$cci) <- "Charlson's comorbidity index"
    levels(df_model$cci) <- c("Low", "Moderate", "Severe")
    var_label(df_model$agegroup) <- "Age group"
    df_model$agegroup <- relevel(df_model$agegroup, ref = "50-59")
    
    # Mediators
    var_label(df_model$bmi2) <- "BMI (centred)"
    var_label(df_model$alc) <- "Alcohol misuse"
    levels(df_model$alc) <- c("No", "Yes")
    var_label(df_model$smokstatus) <- "Smoking status"
    levels(df_model$smokstatus) <- stringr::str_to_title(levels(df_model$smokstatus))
    var_label(df_model$gc90days) <- "Recent high dose glucocorticoid steroid use (<30days)"
    levels(df_model$gc90days) <- c("No", "Yes")
    df_model$sleep <- factor(df_model$sleep, levels = 0:1, labels = c("No", "Yes"))
    var_label(df_model$sleep) <- "Sleep problems"
    var_label(df_model$severity) <- paste(str_to_title(exposure), "severity", sep = " ")
    levels(df_model$severity) <- str_to_title(levels(df_model$severity))
    if (ABBRVexp == "ecz") {
      df_model$comorbid <- df_model$asthma
      var_label(df_model$comorbid) <- "Asthma"
      levels(df_model$comorbid) <- c("No", "Yes")
    } else{
      df_model$comorbid <- df_model$arthritis
      var_label(df_model$comorbid) <- "Arthritis"
      levels(df_model$comorbid) <- c("No", "Yes")
    }
    
    saveRDS(df_model, file = paste0(datapath, "out/nogp/df_model", ABBRVexp, "_", outcome, ".rds"))
  }
}