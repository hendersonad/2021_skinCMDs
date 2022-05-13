library(here)
library(tidyverse)
library(survival)
library(survminer)
library(gt)
library(gtsummary)
library(janitor)
library(flextable)
library(skimr)
library(clipr)
library(KMunicate)
library(data.table)
library(arrow)
library(formattable)
library(labelled)

if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
	  datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

samplingsmall <- F ## set to TRUE if using a small sample of cohort ot build code 
rerunSteroids <- F
XX <- c("psoriasis", "eczema")
export_datasets <- TRUE

for(exposure in XX){
	#exposure <- XX[1]
	ABBRVexp <- str_sub(exposure,1 ,3)
	.dib(exposure) 

	
	# import data -------------------------------------------------------------
	patient <- haven::read_dta(paste0(datapath, "in/Patient_extract_", ABBRVexp, "_extract3_1.dta"))
	cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
	var_depressionDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-definite.dta"))
	var_anxietyDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-definite.dta"))
	var_depressionAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-all.dta"))
	var_anxietyAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-all.dta"))
	var_depressionCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-censor.dta"))
	var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-censor.dta"))
	var_smi <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-smi-definite.dta"))
	var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-age-sex-gp", ".dta"))
	var_bmi <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-BMI-all.dta"))
	var_eth <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-ethnicity.dta"))
	var_alc <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-harmfulalcohol.dta"))
	var_smok <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-smoke-all.dta"))
	var_sev <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-severity.dta"))
	var_sleep_all <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-sleep-all.dta"))
	var_sleep <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-sleep-definite.dta"))
	var_cci <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-cci.dta"))
		# get rid of indexdate in case they get in the way of merging later
		var_bmi <- var_bmi %>% select(-indexdate)
		var_smok <- var_smok %>% select(-indexdate)
		
	death <- patient %>%
		select(patid, tod) 
	cohort_death <- cohort %>%
		left_join(death, by = "patid") 
	cohort_death$death <- cohort_death$tod <= as.numeric(cohort_death$enddate)
	table(cohort_death$death, useNA = "always") ## eczema: 1.4mil NA, 680k died, 134k left for other reasons
	
	if(exposure == "eczema"){
		var_asthma <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-asthma.dta"))
		var_steroids <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-continuousGCRx-90days.dta"))
	}else if(exposure == "psoriasis"){
		var_arthritis <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-rheumatoid-arthritis.dta"))
	}
		
	# get rid of silly stata dbl+lbl format
	var_eth <- var_eth %>% 
		mutate(eth_edited = as_factor(eth5))
			levelsofeth5 <- levels(var_eth$eth_edited) # get all levels of factor var
			levelsofeth5sub <- levelsofeth5[1:6] # get just the levels that start "n. "
			levelsofeth5sub <- str_sub(levelsofeth5sub, 4, str_length(levelsofeth5sub)) # get rid of the "n. "
			levelsofeth5new <- c(levelsofeth5sub, levelsofeth5[7]) # combine
			levels(var_eth$eth_edited) <- levelsofeth5new # assign to data
			var_agesexgp <- var_agesexgp %>% 
			  mutate(gender = as_factor(gender))
			var_smok <- var_smok %>% 
			  mutate(smokstatus = as_factor(smokstatus))
			
			if(exposure == "eczema"){
			  var_sev <- var_sev %>% 
			    mutate(modsevere = as_factor(modsevere))
			}else if(exposure == "psoriasis"){
			  var_sev <- var_sev %>% 
			    rename(date = date_severe) %>%
			    mutate(modsevere = as_factor(src_severetreat))
			  levels(var_sev$modsevere) <- rep("severe",3)
			}
			
			## linked data
			var_catstairs_pt <- readr::read_delim(file = paste0(datapath, "linked/patient_carstairs_20_051.txt"), delim = "\t")
			var_catstairs_prac <- readr::read_delim(file = paste0(datapath, "linked/practice_carstairs_20_051.txt"), delim = "\t")
			var_ruc <- readr::read_delim(file = paste0(datapath, "linked/practice_urban_rural_20_051.txt"), delim = "\t")
			
			# get one variable for RUC
			var_ruc <- var_ruc %>%
			  group_by(pracid) %>%
			  mutate(ruc = max(ni2015_urban_rural, e2011_urban_rural, w2011_urban_rural, s2016_urban_rural, na.rm = T)) %>%
			  ungroup()
			
			# expand steroid risk window data to on/off -------------------------------
			if(exposure=="eczema"){
			  if(rerunSteroids == F){
			    steroids_exp <- readRDS(paste0(datapath, "out/", ABBRVexp, "steroids_exp.rds"))
			  }else{
			    steroids_exp <- var_steroids %>%
			      left_join(select(cohort, patid, indexdate), by = "patid")
			    steroids_exp <- steroids_exp[rep(1:nrow(var_steroids), each = 2), ] #Duplicate each row
			    steroids_exp[1:nrow(steroids_exp) %% 2 == 0, ] <- NA #Set every second row to NA
			    steroids_exp <- fill(steroids_exp, indexdate, patid, .direction = "down") #Copy values down to empty cells
			    steroids_exp$active <- rep(c("active", "not active"), nrow(var_steroids)) #Make every first row "active", and every second row "not active"
			    
			    #Make starts and ends for the non active periods
			    #test <- steroids_exp %>% slice(1:100)
			    
			    steroids_exp <- steroids_exp %>% 
			      group_by(patid) %>%
			      mutate(end=case_when(!is.na(enddate) ~ as.Date(enddate), #If its already there
			                           is.na(lead(startdate)) ~ as.Date("3000-01-01"), #If its a participants last observation
			                           is.na(enddate) ~ as.Date(lead(startdate))), #Else, make it the start date of the next observation
			             start=case_when(!is.na(startdate) ~ as.Date(startdate),
			                             is.na(startdate) ~ as.Date(lag(enddate)))) %>%
			      filter(start<end & end>=indexdate) #Keep only observations where the start is before the end and that end within or after the follow up period
			    
			    saveRDS(steroids_exp, file = paste0(datapath, "out/", ABBRVexp, "steroids_exp.rds"))
			  }
			}
			
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
			
			# deal with duplicate patids ----------------------------------------------
			## people who contribute time as a control and a case 
			# cohort_obs <- cohort %>% 
			# 	group_by(patid) %>% 
			# 	mutate(obs = 1:n()) 
			# ## get those that are not duplicated
			# cohort_nodups <- cohort_obs %>%
			# 	filter(obs==1)
			# ## get duplciates and rename patid for exposed period 
			# cohort_dups <- cohort_obs %>%
			# 	filter(max(obs)>1) %>%
			# 	ungroup() %>%
			# 	filter(exposed==0) %>%
			# 	mutate(newpatid = as.numeric(paste0("9999", patid))) %>%
			# 	select(-patid) %>%
			# 	select(patid = newpatid, everything())
			# 
			# cohort_edited <- cohort_nodups %>%
			# 	bind_rows(cohort_dups)
			# if(nrow(cohort) != nrow(cohort_edited)){
			# 	.dib("You've messed up")
			# }else{
			# 	.dib("New df same size as original...carry on")
			# }
			# x <- sum(duplicated(cohort_edited$patid))
			# if(x!=0){
			# 	.dib("Still got some duplication in your data pal")
			# }else{
			# 	.dib(paste0("Got rid of ", prettyNum(sum(duplicated(cohort$patid))), " duplicates"))
			# }
			
			# make a decent sample of the main cohort ---------------------------------
			if(samplingsmall == TRUE){
			  set.seed(105351)
			  sample_sets <- sample(cohort$setid, size = 1000)
			  
			  cohort_sample <- cohort %>%
			    filter(setid %in% sample_sets)
			  ## Store the edited cohort sample 
			  if(export_datasets){
			    saveRDS(cohort_sample, file = paste0(datapath, "out/", ABBRVexp, "_cohort_sample.rds"))
			  }
			  
			  cohort_edited <- cohort_sample
			  cohort_edited <- cohort_edited %>%
			    select(setid, patid, exposed, indexdate, enddate) %>%
			    left_join(var_agesexgp, by = "patid")
			  export_datasets <- FALSE
			}else{
			  cohort_edited <- cohort
			  cohort_edited <- cohort_edited %>%
			    select(setid, patid, exposed, indexdate, enddate) %>%
			    left_join(var_agesexgp, by = "patid")
			}
			
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
			          file = paste0(datapath, "out/", ABBRVexp, "_an_anxiety.rds"))
			}
			if(sum(grepl(pattern = "cohort_edited", x = ls()))==0){
			  an_anxiety <- readRDS(paste0(datapath, "out/", ABBRVexp, "_an_anxiety.rds"))
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
			          file = paste0(datapath, "out/", ABBRVexp, "_an_depression.rds"))
			}
			
			if(sum(grepl(pattern = "cohort_edited", x = ls()))==0){
			  an_depression <- readRDS(paste0(datapath, "out/", ABBRVexp, "_an_depression.rds"))
			}
			
			
			# tmerge function ---------------------------------------------------------
			build_tmerge <- function(exp = 1, df_outcome = an_anxiety, outcome = var_anx){
			  df_exp <- df_outcome %>%
			    filter(exposed==exp)
			  
			  out_tmerge0 <- tmerge(df_exp, df_exp, id=patid, tstart = indexdate, tstop = enddate) #This first step sets the range of follow up
			  out_tmerge1 <- tmerge(out_tmerge0, outcome, id=patid, out = event(eventdate))
			  out_tmerge2 <- tmerge(out_tmerge1, var_alc, id=patid, alc = tdc(eventdate)) 
			  out_tmerge3 <- tmerge(out_tmerge2, var_sleep, id=patid, sleep = tdc(eventdate)) 
			  out_tmerge3b <- tmerge(out_tmerge3, var_smok, id=patid, smokstatus = tdc(eventdate, smokstatus, exp)) 
			  if(exposure == "eczema"){
			    out_tmerge4 <- tmerge(out_tmerge3b, var_sev, id=patid, severity = tdc(date, modsevere, exp))
			    out_tmerge5 <- tmerge(out_tmerge4, var_asthma, id=patid, asthma = tdc(eventdate)) 
			    out_tmerge6 <- tmerge(out_tmerge5, steroids_exp, id=patid, gc90days = tdc(start, active, "not active"))
			  }else if(exposure == "psoriasis"){
			    out_tmerge4 <- tmerge(out_tmerge3b, var_sev, id=patid, severity = tdc(date, modsevere, exp))
			    out_tmerge5 <- tmerge(out_tmerge4, var_arthritis, id=patid, arthritis = tdc(eventdate)) 
			    out_tmerge6 <- out_tmerge5 ## don't have steroids but it will break the rest of the function code if there's no out_tmerge6
			  }
			  # all sleep
			  out_tmerge7 <- tmerge(out_tmerge6, var_sleep_all, id=patid, sleep_all = tdc(eventdate)) 
			  # censoring events
			  out_tmerge8 <- tmerge(out_tmerge7, death, id=patid, death = event(tod))
			  #out_tmerge9 <- tmerge(out_tmerge8, var_smi, id=patid, death = tdc(eventdate)) ## not going to work because of manual censoring above 
			  
			  ## Restrict to first occurence of outcome ## post-processing to keep only records upto event if have the event. 
			  out_tmerge9 <- out_tmerge8 %>%
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
			  saveRDS(anxiety_full, file = paste0(datapath, "out/", ABBRVexp, "-anxiety_full.rds")) ## some temp stats from wrong run of tmerge commands - 3541262 (16 cols). Now 7234333 (because I forgot to change exp=0 in the unexposed run like a fool)
			  saveRDS(depression_full, file = paste0(datapath, "out/", ABBRVexp, "-depression_full.rds")) ## some temp stats from wrong run of tmerge commands - 6286583 (16 cols)		
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
			rm(an_anxiety, an_depression, anxiety_exposed, anxiety_split, anxiety_unexposed, depression_exposed, depression_split, depression_unexposed,
			   cohort, cohort_edited, 
			   KM)
			
			## temp fix - delete later when run full code
			# depression_full <- depression_full %>%
			# 	mutate(gender=as_factor(gender)) %>% 
			# 	mutate(severity = as_factor(severity))
			# anxiety_full <- anxiety_full %>%
			# 	mutate(gender=as_factor(gender)) %>% 
			# 	mutate(severity = as_factor(severity))
			
			if(sum(grepl(pattern = "anxiety_full", x = ls()))==0){
			  anxiety_full <- readRDS(paste0(datapath, "out/", ABBRVexp, "-anxiety_full.rds"))
			  depression_full <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_full.rds"))
			  anxiety_full$dob <- as.Date(paste(anxiety_full$realyob, "06-01", sep="-"))
			  depression_full$dob <- as.Date(paste(depression_full$realyob, "06-01", sep="-"))
			}
			
			merge_static <- function(df_out = depression_full){
			  if(exposure == "eczema"){
			    df_static <- df_out %>% 
			      mutate(arthritis = NA) %>%
			      select(setid, patid, pracid, exposed, indexdate, enddate, dob, gender, tstart, tstop, out, asthma, arthritis, alc, smokstatus, severity, sleep, sleep_all, gc90days, death) 
			  }else if(exposure == "psoriasis"){
			    df_static <- df_out %>% 
			      mutate(gc90days = NA, asthma = NA) %>%
			      select(setid, patid, pracid, exposed, indexdate, enddate, dob, gender, tstart, tstop, out, asthma, arthritis, alc, smokstatus, severity, sleep, sleep_all, gc90days, death) 
			  }
			  df_static <- df_static %>% 
			    left_join(select(var_eth, patid, eth_edited), by="patid") %>%
			    left_join(var_bmi, by="patid") %>%
			    #left_join(var_smok, by="patid") %>%
			    #rename(smokdate = eventdate) %>%
			    left_join(select(var_catstairs_pt, patid, carstairsPt=carstairs2011_5), by="patid") %>%
			    left_join(select(var_catstairs_prac, carstairsPrac=carstairs2011_5, everything()), by="pracid") %>%
			    left_join(select(var_ruc, pracid, ruc), by="pracid") %>%
			    mutate(carstairs = carstairsPt) %>%
			    mutate_at("carstairs", ~ifelse(is.na(.), carstairsPrac, .)) %>%
			    select(-carstairsPrac, -carstairsPt) %>% 
			    left_join(select(var_cci, patid, cci), by="patid") 
			  
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
			                                         "66+")),
			      bmi_cat=case_when(bmi < 18.5 ~ "underweight",
			                        bmi < 25 ~ "normal",
			                        bmi < 30 ~ "overweight",
			                        bmi >= 30 ~ "obese"),
			      bmi_cat=factor(bmi_cat, levels = c("underweight", 
			                                         "normal", 
			                                         "overweight", 
			                                         "obese")),
			      alc=factor(alc),
			      asthma=factor(asthma),
			      arthritis=factor(arthritis),
			      ruc=factor(ruc),
			      carstairs=factor(carstairs),
			      cci=factor(cci)
			    )
			  levels(df_static$ruc) <- c("urban", "rural")
			  levels(df_static$smokstatus)[5:6] <- NA
			  
			  if(exposure == "eczema"){
			    levels(df_static$severity) <- c("moderate", "severe", "mild", "none")
			    df_static$severity <- factor(df_static$severity, levels = c("none", "mild", "moderate", "severe"))
			  }else if(exposure == "psoriasis"){
			    levels(df_static$severity) <- c("severe", "mild", "none")
			    df_static$severity <- factor(df_static$severity, levels = c("none", "mild", "severe"))
			  }
			  
			  df_static
			}
			anxiety_static <- merge_static(anxiety_full)
			depression_static <- merge_static(depression_full)
			
			if(export_datasets){
			  saveRDS(anxiety_static, file = paste0(datapath, "out/", ABBRVexp, "-anxiety_static.rds")) 
			  saveRDS(depression_static, file = paste0(datapath, "out/", ABBRVexp, "-depression_static.rds")) 
			}
			if(sum(grepl(pattern = "anxiety_static", x = ls()))==0){
			  anxiety_static <- readRDS(paste0(datapath, "out/", ABBRVexp, "-anxiety_static.rds"))
			  depression_static <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_static.rds"))
			}
			
			build_split <- function(df_in = anxiety_static){
			  # tosplit <- df_in %>%
			  # 	arrange(setid, patid) %>%
			  # 	slice(1:2000)
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
			rm(anxiety_full, depression_full, var_agesexgp, var_eth, var_bmi, var_smok, var_ruc, var_catstairs_prac, var_catstairs_pt, var_cci)
			
			anxiety_split <- build_split(anxiety_static)
			depression_split <- build_split(depression_static)
			
			table(anxiety_split$severity)
			if(export_datasets){
			  saveRDS(anxiety_split, file = paste0(datapath, "out/", ABBRVexp, "-anxiety_split.rds")) 
			  saveRDS(depression_split, file = paste0(datapath, "out/", ABBRVexp, "-depression_split.rds")) 
			}
			
			
			samplingsmall <- F ## set to TRUE if using a small sample of cohort ot build code 
			rerunSteroids <- F
			XX <- c("psoriasis", "eczema")
			export_datasets <- TRUE
}
for(exposure in XX){
  #exposure <- XX[1]
  ABBRVexp <- str_sub(exposure,1 ,3)
  .dib(exposure) 
  
  anxiety_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-anxiety_split.rds"))
  depression_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_split.rds"))
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
		
		## impute smoking data 
		df_model$smokstatus %>% table(useNA = "always")
		df_model <- df_model %>% 
		  group_by(setid, patid) %>% 
		  mutate(smok_missing_always = ifelse(any(!is.na(smokstatus)), 0, 1),
		         bmi_missing_always = ifelse(any(!is.na(bmi2)), 0, 1)
		  )
		df_model$smokstatus_nomiss <- df_model$smokstatus
		df_model$smokstatus_nomiss[is.na(df_model$smokstatus_nomiss) & df_model$smok_missing_always == 0] <- "Non-Smoker"
		df_model$smokstatus_nomiss %>% table(useNA = "always")
		
		df_model$bmi_miss <- 0
		df_model$bmi_miss[is.na(df_model$bmi2)] <- 1
		df_model$bmi_miss %>% table(useNA = "always")
		df_model$bmi2[is.na(df_model$bmi2) & df_model$bmi_missing_always == 0] <- 0
		
		# re-categorise BMI
		df_model <- df_model %>% 
		  ungroup() %>% 
		  mutate(obese_cat=case_when(bmi < 30 ~ "no evidence of obesity",
		                           bmi < 35 ~ "obese class I (30-34.9)",
		                           bmi < 40 ~ "obese class II (35-39.9)",
		                           bmi >= 40 ~ "obese class III (40+)"))
		df_model$obese_cat <- factor(df_model$obese_cat, levels = c("no evidence of obesity", 
		                                            "obese class I (30-34.9)", 
		                                            "obese class II (35-39.9)", 
		                                            "obese class III (40+)"))
		df_model$obese_cat %>% table(useNA = "always")
		df_model$obese_cat[is.na(df_model$bmi)] <- "no evidence of obesity"
		df_model$obese_cat %>% table(useNA = "always")
		
		saveRDS(df_model, file = paste0(datapath, "out/df_model", ABBRVexp, "_", outcome, ".rds"))
	}
}

