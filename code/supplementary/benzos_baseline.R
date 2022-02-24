library(tidyverse)
library(here)
library(janitor)
library(timetk)
library(gtsummary)
library(magrittr)
library(survival)
library(survminer)
library(htmltools)
library(lubridate)

if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
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

samplingSmall <- F
if(samplingSmall){
	df_anx_split <- df_anx_split %>% slice(1:1e6)
	df_dep_split <- df_dep_split %>% slice(1:1e6)
}
df_anx_split$t <- as.numeric(df_anx_split$tstop - df_anx_split$tstart)
df_dep_split$t <- as.numeric(df_dep_split$tstop - df_dep_split$tstart)

df_anx_split$gc90days <- factor(df_anx_split$gc90days, levels = c("not active", "active"))
df_dep_split$gc90days <- factor(df_dep_split$gc90days, levels = c("not active", "active"))

df_anx_split$bmi2 <- df_anx_split$bmi - mean(df_anx_split$bmi, na.rm = T)
df_dep_split$bmi2 <- df_dep_split$bmi - mean(df_dep_split$bmi, na.rm = T)



var_sleep_all <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-sleep-all.dta"))
var_sleep <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-sleep-definite.dta"))
patient <- haven::read_dta(paste0(datapath, "in/Patient_extract_", ABBRVexp, "_extract3_1.dta"))
cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-age-sex-gp", ".dta"))
var_depressionDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-definite.dta"))
var_anxietyDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-definite.dta"))
# might look at deprivation too but going to be missing a fair bit
var_catstairs_pt <- readr::read_delim(file = paste0(datapath, "linked/patient_carstairs_20_051.txt"), delim = "\t")
var_catstairs_prac <- readr::read_delim(file = paste0(datapath, "linked/practice_carstairs_20_051.txt"), delim = "\t")

# some headline stats -----------------------------------------------------
df_all <- patient %>%
	select(patid, gender, realyob) %>%
	left_join(var_sleep_all, by = "patid") %>%
	rename(sleepdate = eventdate) %>%
	select(patid, gender, realyob, sleepdate, ingredient, bnftext, src) %>%
	left_join(select(cohort, setid, patid, exposed, indexdate, enddate)) %>%
	left_join(select(var_anxietyDef, patid, eventdate), by = "patid") %>%
	rename(anxiety = eventdate) %>%
	left_join(select(var_depressionDef, patid, eventdate), by = "patid") %>%
	rename(depression = eventdate) %>%
	mutate(
		sleep = case_when(
			!is.na(sleepdate) ~ 1,
			is.na(sleepdate) ~ 0),
		anx = case_when(
			!is.na(anxiety) ~ 1,
			is.na(anxiety) ~ 0),
		dep = case_when(
			!is.na(depression) ~ 1,
			is.na(depression) ~ 0),
		) %>%
	mutate_at(c("anx", "dep", "sleep"), ~as.factor(.))

skimr::skim(df_all)
# 6,022,907 participants, 1.6m with depression. 1.18m with anxiety, 1.85m with sleep

df_all %>%
	tabyl(anx, dep) %>%
	adorn_percentages("all") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## so 13% of people have both anx and dep. 6% anx only and 14% depression only

## by sleep status
df_all %>%
	tabyl(anx, dep, sleep) %>%
	adorn_percentages("all") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns()	## among people with a sleep prescription: 28% both, 10% anxiety only, 21% depression only
							## among people WITHOUT a sleep prescription: 6.3% both, 5% anxiety only, 11% depression only
# So, are sleep meds given more often in people with anxiety (to treat sleep and mental illness perhaps)? Probably not

# benzo prescriptions -----------------------------------------------------
## who is getting benzos and definite sleep meds
df_sleep <- df_all %>% 
	left_join(select(var_sleep, patid, sleep_defdate = eventdate, ingredient_def = ingredient, bnftextdef = bnftext, src_def = src),
						by = "patid") %>%
	mutate(sleep_def = as.factor(case_when(
		!is.na(sleep_defdate) ~ 1,
		is.na(sleep_defdate) ~ 0
	)))
df_sleep %>% 
	tabyl(sleep_def, sleep) %>% 
	adorn_percentages("row") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## so 16% of people without Definite sleep code do have a benzo code

df_sleep %>% 
	tabyl(sleep_def, sleep) %>% 
	adorn_percentages("all") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## and 18% (1mil) have both sleep codes, another 13% have benzo codes. Only 69% of pop do not have sleep code

df_ingredient <- df_sleep %>% 
	filter(sleep == 1) %>%
	#select(patid, exposed, realyob, indexdate, sleep, sleep_def, anx, dep, ingredient, ingredient_def, src, src_def)
	mutate(dob = as.Date(paste(realyob, "07", "01", sep = "-"))) 
df_ingredient <- df_ingredient %>% 
	mutate(age_at_sleep = as.numeric(sleepdate - dob)/365.25) %>%
	mutate(exposed = as.factor(exposed))

df_ingredient %>% 
	tabyl(src_def, src) %>%
	adorn_totals(c("row", "col")) %>%
	adorn_percentages("row") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## 19% of sleep codes come from read (all definites), 81% from prescriptions

df_ingredient %>%
	tabyl(ingredient, exposed) %>% 
	adorn_totals("col") %>%
	adorn_percentages("col") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## 37% of control codes are diazepam. Only 34% in eczema

ggplot(df_ingredient, aes(x = age_at_sleep, group = exposed, col = exposed, fill = exposed)) +
	geom_histogram(alpha = 0.2) +
	facet_wrap(~sleep_def) +
	theme_classic()


# filter on baseline ------------------------------------------------------
df_sleep %>%
	mutate(sleep_censor = as.numeric(sleep == 1 & sleepdate <= indexdate)) %>%
	tabyl(exposed, sleep_censor) %>% 
	adorn_totals("row") %>%
	adorn_percentages("row") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## at baseline. 26% of eczema have sleep codes and 17% of controls

df_sleep %>%
	mutate(sleep_censor = as.numeric(sleep == 1 & sleepdate <= indexdate),
				 sleep_benzo = as.numeric(sleep_censor == 1 & sleep_def == 0)) %>%
	tabyl(exposed, sleep_benzo) %>% 
	adorn_totals("row") %>%
	adorn_percentages("row") %>%
	adorn_pct_formatting(digits = 2) %>%
	adorn_ns() ## at baseline. 9% of ecemza have benzo and 7% of controls
