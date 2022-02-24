pacman::p_load(here, haven, tidyverse, ggplot2, survival, survminer)
pacman::p_load(gt, gtsummary)
pacman::p_load(janitor, flextable)
#pacman::p_load(summarytools)
pacman::p_load(skimr)
pacman::p_load(clipr)
pacman::p_load(KMunicate)
pacman::p_load(data.table)
pacman::p_load(arrow)

datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"

exposure <- "eczema"
ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure) 

# import data -------------------------------------------------------------
cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-", exposure, "-main-mhealth.dta"))
patinfo <- haven::read_dta(paste0(datapath, "in/Patient_extract_", ABBRVexp, "_extract3_1.dta"))
var_depressionDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-definite.dta"))
var_anxietyDef <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-definite.dta"))
var_depressionAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-all.dta"))
var_anxietyAll <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-all.dta"))
var_depressionCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-depression-censor.dta"))
var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-anxiety-censor.dta"))
var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-age-sex-gp", ".dta"))
var_bmi <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-BMI-all.dta"))
var_eth <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-ethnicity.dta"))
var_alc <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-harmfulalcohol.dta"))
var_smok <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-smoke-all.dta"))
var_sev <- haven::read_dta(paste0(datapath, "out/variables-", ABBRVexp, "-severity.dta"))

# get rid of indexdate because they get in the way of merging later
var_bmi <- var_bmi %>% select(-indexdate)
var_smok <- var_smok %>% select(-indexdate)

# merge into one big ol’ dataset ------------------------------------------
## depression 
df_tab1 <- cohort %>%
	left_join(var_depressionCensor, by = "patid") %>%
	#select(-constype, -readcode, -readterm) %>%
	rename(censorDepDate = eventdate) %>%
	mutate_at("censorDepDate", ~ifelse(is.na(.), as.Date("3000-01-01"), .)) %>%
	left_join(var_depressionAll, by = "patid") %>%
	#select(-constype, -readcode, -readterm) %>%
	mutate(dep = case_when(
		eventdate < censorDepDate &
		eventdate < enddate &
		eventdate > indexdate &
		!is.na(eventdate) ~ 1,
		TRUE ~ 0
	)) %>%
	rename(depdate = eventdate) 

df_tab1 %>% 
	filter(!is.na(depdate)) %>%
	filter(depdate<censorDepDate) %>%
	filter(depdate<enddate) %>%
	filter(depdate>=indexdate) %>%
	count()

## anxiety 
df_tab1 <- df_tab1 %>%
	left_join(var_anxietyCensor, by = "patid") %>%
	#select(-constype, -readcode, -readterm) %>%
	rename(censorAnxDate = eventdate) %>%
	mutate_at("censorAnxDate", ~ifelse(is.na(.), as.Date("3000-01-01"), .)) %>%
	left_join(var_anxietyAll, by = "patid") %>%
	#select(-constype, -readcode, -readterm) %>%
	mutate(anx = case_when(
		eventdate < censorAnxDate & 
		eventdate < enddate &
		eventdate > indexdate &
		!is.na(eventdate) ~ 1,
		TRUE ~ 0
	)) %>%
	rename(anxdate = eventdate)

## demographics and lifestyle vars
df_tab1 <- df_tab1 %>%
	left_join(var_agesexgp, by = "patid") %>%
	left_join(var_eth, by="patid") %>%
	left_join(var_bmi, by="patid") %>%
	mutate_at("bmi", ~ifelse(dobmi>enddate, NA, .)) %>%
	left_join(var_alc, by="patid") %>%
	#select(-ingredient, -prodcode, -bnftext, -src, -medcode, -readterm) %>%
	rename(alcdate = eventdate) %>%
	mutate_at("alcdate", ~ifelse(.>enddate, NA, .)) %>%
	mutate(alc = ifelse(is.na(alcdate), 0, 1)) %>%
	left_join(var_smok, by="patid") %>%
	rename(smokdate = eventdate) %>%
	mutate_at("smokstatus", ~ifelse(smokdate>enddate, NA, .))

# table(is.na(as.vector(df_tab1$bmi)))
# table(df_tab1$alc)
# table(df_tab1$smokstatus)
# sum(is.na(df_tab1$smokstatus))

## add in severity data (but this will depend on exposure)
if(str_detect(str_to_lower(exposure),"eczema")){
	df_tab1 <- df_tab1 %>%
	left_join(var_sev, by = "patid") %>%
	mutate_at("modsevere" , ~replace_na(., 0)) %>%
	rename(date_severe = date) %>%
	rename(severity = modsevere)
}
if(str_detect(str_to_lower(exposure),"psoriasis")){
	df_tab1 <- df_tab1 %>%
	left_join(var_sev, by = "patid") %>%
	mutate_at("src_severetreat" , ~ifelse(date_severe < indexdate, 0, .)) %>%
	mutate_at("src_severetreat" , ~ifelse(date_severe > enddate, 0, .)) %>%
	mutate(severe = ifelse(is.na(date_severe), 0, 1)) %>%
	rename(severity = src_severetreat)
}

glimpse(df_tab1)

write_parquet(df_tab1, sink = paste0(datapath, "out/df_tab1_",exposure,".gz.parquet"))
#df_tab1 <- read_parquet(file = paste0(datapath, "out/df_tab1_",exposure,".gz.parquet"))

df_tab1 %>%
	count(dep)
df_tab1 %>%
	count(anx)

## get rid of big read in data to save memory
rm(list=c(
	"patinfo",
	"var_depressionDef",
	"var_anxietyDef",
	"var_agesexgp",
	"var_bmi",
	"var_eth",
	"var_alc",
	"var_smok",
	"var_sev"))

# bit of formatting 
df_tab1 <- df_tab1 %>%
	mutate(dob = as.Date(paste(realyob, "01", "01", sep = "-")),
				 age = as.numeric((indexdate - dob)/365.25)) %>%
	mutate_at(c("exposed", "gender", "dep", "anx", "eth5", "eth16", "alc", "smokstatus", "severity"), ~as_factor(.))

## Generate new variables
# follow up time 
# age cat
df_tab1 <- df_tab1 %>%
	mutate(futime = as.numeric((enddate-indexdate)/365.25),
				 age_cat = case_when(
				 	age < 30 ~ "18 - 29",
				 	age < 40 ~ "30 - 39",
				 	age < 50 ~ "40 - 49",
				 	age < 66 ~ "50 - 65",
				 	age >= 66 ~ "66+"),
				 bmi_cat = case_when(
				 	bmi < 20 ~ "<20",
				 	bmi < 25 ~ "20-24",
				 	bmi < 30 ~ "25-29",
				 	bmi < 66 ~ "30+")
	) %>%
	mutate_at(c("age_cat","bmi_cat"), ~as_factor(.))

# summarise using different options ---------------------------------------
## base R
#summary(df_tab1)

## skimr
#skim(df_tab1)
#group_by(df_tab1, exposed) %>% skim()

## summarytools
#summarytools::descr(df_tab1)
#summarytools::dfSummary(df_tab1)

tab1 <- df_tab1 %>% 
	select(patid, exposed, futime, gender, age, age_cat, eth5, bmi, bmi_cat,
				 dep, anx, alc, smokstatus, severity) %>%
	mutate_at("exposed", ~ifelse(. == 0, "Matched controls", paste0("With ", exposure)))

## reorder value labels 
tab1 <- tab1 %>% 
	mutate(gender = factor(gender, levels = c("Male", "Female"), exclude = c("Data Not Entered", "Indeterminate", "Unknown")),
				 age_cat = factor(age_cat, levels = c("18 - 29", "30 - 39", "40 - 49", "50 - 65", "66+")),
				 bmi_cat = factor(bmi_cat, levels = c("<20", "20-24", "25-29", "30+")),
				 eth5 = factor(eth5, 
				 							levels = c("0. White", "1. South Asian", "2. Black", "3. Other","4. Mixed", "5. Not Stated"),
				 							labels = c("White", "South Asian", "Black", "Other", "Mixed", "Not Stated"), exclude = "equally common"),
				 alc = factor(alc, levels = 0:1, labels = c("No", "Yes")),
				 dep = factor(dep, levels = 0:1, labels = c("No", "Yes")),
				 anx = factor(anx, levels = 0:1, labels = c("No", "Yes")))

## add in severity data (but this will depend on exposure)
if(str_detect(str_to_lower(exposure),"eczema")){
	tab1 <- tab1 %>%
		mutate(severity = factor(severity, levels = levels(tab1$severity), labels = c("Mild", "Moderate", "Severe")))
	
}
if(str_detect(str_to_lower(exposure),"psoriasis")){
	tab1 <- tab1 %>%
		mutate(severity = factor(severity, levels = 0:3, labels = c("None", "Systemic", "Phototherapy", "Biologic")))
}

## investigating the unexposed people with severe psoriasis Rx codes
table(tab1$exposed, tab1$severity)

unexposed_severe <- tab1 %>%
	filter(severity !="None" & exposed == "Matched controls") %>%
	pull(patid)
unexposed_severe <- as.vector(unexposed_severe)

df_unexp_severe <- df_tab1 %>% 
	filter(patid %in% unexposed_severe) %>%
	select(setid, patid, exposed, enddate, indexdate) %>%
	group_by(patid) %>%
	mutate(count = 1:n(),
				 obs = max(count)) %>%
	ungroup() 
df_unexp_severe %>%
	count(obs) ## so almost all are poeple who contribute some unexposed time as a control, then become exposed
df_tab1 %>% 
	filter(patid == 269504)
unexp_problems <- df_unexp_severe %>%
	filter(obs == 1) %>%
	pull(patid) %>%
	as.vector()
df_problems <- df_tab1 %>%
	filter(patid %in% unexp_problems) ## but there are 4 records of people who have  severe Rx but never met the full Psoriasis exposed definition

## remove some stuff 
rm(df_tab1)

## make table with nice gtsummary 
table1 <- tab1 %>% 
	select(-patid) %>%
	tbl_summary(by = exposed,
							statistic = list(all_continuous() ~ "{p50} ({p25}-{p75})",
															 #futime ~ "{sum} ({p50}; {p25}-{p75})",
															 all_categorical() ~ "{n} ({p}%)"),
							digits = all_continuous() ~ 1,
							label = list(gender = "Sex",
													 futime = "Follow-up time (years)",
													 age = "Age",
													 age_cat = "Age (categorised)",
													 eth5 = "Ethnicity", 
													 bmi = "BMI",
													 bmi_cat = "BMI (categorised)",
													 dep = "Depression",
													 anx = "Anxiety", 
													 alc = "Harmful alcohol use", 
													 smokstatus = "Smoking status",
													 severe = "Severe")
	) %>%
	add_overall() %>%
	bold_labels() %>%
	modify_footnote(
		all_stat_cols() ~ "Median (IQR) or Frequency (%)"
	) 

table1

table1 %>%
	as_flex_table() %>%
	flextable::save_as_docx(path = here::here("out", paste0("tab1_",ABBRVexp,".docx")))
