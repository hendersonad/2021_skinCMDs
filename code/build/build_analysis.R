pacman::p_load(haven, tidyverse, ggplot2, survival, survminer)

datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2020_eczema_extract/"

cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-eczema-main-mhealth.dta"))
var_depressionDef <- haven::read_dta(paste0(datapath, "out/ecz-depression-definite.dta"))
var_anxietyDef <- haven::read_dta(paste0(datapath, "out/ecz-anxiety-definite.dta"))
## All 
var_depressionAll <- haven::read_dta(paste0(datapath, "out/ecz-depression-all.dta"))
var_anxietyAll <- haven::read_dta(paste0(datapath, "out/ecz-anxiety-all.dta"))

## Static vars
var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-ecz-age-sex-gp", ".dta"))
var_bmi <- haven::read_dta(paste0(datapath, "out/ecz-BMI-all.dta"))
var_eth <- haven::read_dta(paste0(datapath, "out/ecz-ethnicity.dta"))
var_alc <- haven::read_dta(paste0(datapath, "out/ecz-harmfulalcohol.dta"))
var_smok <- haven::read_dta(paste0(datapath, "out/ecz-smoke-all.dta"))

# get rid of indexdate because they get in the way of merging later
var_bmi <- var_bmi %>% select(-indexdate)
var_smok <- var_smok %>% select(-indexdate)

# get a small sample to make life a little easier  ------------------------
cases <- cohort$patid[cohort$exposed == 1]
set.seed(0153825)
#case_sample <- sample(size = 1000, x = cases) %>%
#	as_tibble()
case_sample <- cohort$patid %>%
	as_tibble()

cohort_sample <- cohort %>%
	inner_join(case_sample, by = c("setid" = "value"))

# get their DEFINITE outcome date (anxiety or depression) -----------------
add_in_var <- function(df, df_var, shorthand = "dep"){
	shorthand <- as.character(shorthand)
	newvar <- paste0(shorthand, "date")
	cohort_sample_outcome <- df %>%
		left_join(df_var, by = "patid") %>%
		mutate(binary = ifelse(is.na(eventdate), 0, 1))
	names_vec <- names(cohort_sample_outcome)
	names_vec[names_vec=="binary"] <- shorthand
	names_vec[names_vec=="eventdate"] <- newvar
	names(cohort_sample_outcome) <- names_vec
	return(cohort_sample_outcome)
}
cohort_sample_dep <- add_in_var(df = cohort_sample, df_var = var_depressionDef, shorthand = "dep")
cohort_sample_depanx <- add_in_var(df = cohort_sample_dep, df_var = var_anxietyDef, shorthand = "anx")
table(cohort_sample_depanx$dep)
table(cohort_sample_depanx$anx)

# censor on outcome (new enddate) for schiz and BPD as well ---------------
cohort_sample_censor <- cohort_sample_depanx %>%
	mutate(censordate = as.Date("3000-01-01")) %>%
	mutate_at("censordate", ~pmin(censordate, depdate, anxdate , na.rm = T)) %>%
	mutate(newenddate = pmin(enddate, censordate)) %>%
	dplyr::select(-censordate, -enddate, enddate = newenddate) %>%
	mutate(censor = ifelse(enddate<indexdate, 1, 0)) 

cohort_sample_censor <- cohort_sample_censor %>%
	filter(exposed == 1 | exposed == 0 & censor == 0)

# Get rid of invalid sets (no case)  --------------------------------------
caseid_censor <- cohort_sample_censor %>%
	filter(exposed == 1 & censor == 1) %>%
	select(patid, casecensor = censor)

cohort_sample_valid <- cohort_sample_censor %>%
	left_join(caseid_censor, by = c("setid" = "patid")) %>%
	filter(is.na(casecensor))

length(unique(cohort_sample_depanx$setid)) ## 1000
length(unique(cohort_sample_valid$setid)) ## 749

# merge in outcomes -------------------------------------------------------
analysis_sample <- cohort_sample_valid %>%
	select(-depdate, -dep, -anxdate, -anx, -casecensor, -censor)
analysis_sample_dep <- add_in_var(df = analysis_sample, df_var = var_depressionAll, shorthand = "dep")
analysis_sample_depanx <- add_in_var(df = analysis_sample_dep, df_var = var_anxietyAll, shorthand = "anx")

table(analysis_sample_dep$dep)
table(analysis_sample_depanx$anx)

## time updated vars 

## add in static vars - BMI, alc, ethnicity, smoking, dob and sex
analysis_dta <- analysis_sample_depanx

analysis_dta <- analysis_dta %>% 
	left_join(var_agesexgp, by="patid") 
analysis_dta <- analysis_dta %>% 
	left_join(var_bmi, by="patid") 
analysis_dta <- analysis_dta %>% 
	left_join(var_eth, by="patid") 
analysis_dta <- analysis_dta %>% 
	left_join(var_alc, by="patid") 
analysis_dta <- analysis_dta %>% 
	left_join(var_smok, by="patid") 


# run an easy model  ------------------------------------------------------
analysis_dta <- analysis_dta %>%
	mutate(time = enddate-indexdate)
anx_km <- survfit(Surv(time,anx)~as.factor(exposed), data=analysis_dta)
plot(anx_km,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,1,c("None","Exposed"),col=c("black","grey"),lty=1,cex=0.5)

with(analysis_dta, table(exposed, anx))

dep_km <- survfit(Surv(time,dep)~as.factor(exposed), data=analysis_dta)
plot(dep_km,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,1,c("None","Exposed"),col=c("black","grey"),lty=1,cex=0.5)

with(analysis_dta, table(exposed, dep))

dep_cox<-coxph(Surv(time,dep)~as.factor(exposed),data=analysis_dta)
summary(dep_cox)

dep_survfit <- survfit(dep_cox,newdata=data.frame(exposed=c(0,1)))
summary(dep_survfit)

plot(dep_survfit,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,1,c("Non","Eczema"),col=c("black","grey"),lty=1,cex=0.5)

## add in comorbidities? 

# Identify earliest index and latest enddate ------------------------------


# stset the data ----------------------------------------------------------


# Split on age ------------------------------------------------------------


# Split on calendar time --------------------------------------------------


# Save --------------------------------------------------------------------





