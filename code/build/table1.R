pacman::p_load(here, haven, tidyverse, ggplot2, survival, survminer)
pacman::p_load(gt, gtsummary)
pacman::p_load(janitor, flextable)
#pacman::p_load(summarytools)
pacman::p_load(skimr)
pacman::p_load(clipr)
pacman::p_load(KMunicate)

datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2020_eczema_extract/"

cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-eczema-main-mhealth.dta"))
patinfo <- haven::read_dta(paste0(datapath, "in/Patient_extract_ecz_extract3_1.dta"))
var_depressionDef <- haven::read_dta(paste0(datapath, "out/ecz-depression-definite.dta"))
var_anxietyDef <- haven::read_dta(paste0(datapath, "out/ecz-anxiety-definite.dta"))
var_depressionCensor <- haven::read_dta(paste0(datapath, "out/ecz-depression-censor.dta"))
var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/ecz-anxiety-censor.dta"))
var_depressionAll <- haven::read_dta(paste0(datapath, "out/ecz-depression-all.dta"))
var_anxietyAll <- haven::read_dta(paste0(datapath, "out/ecz-anxiety-all.dta"))
var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-ecz-age-sex-gp", ".dta"))
var_bmi <- haven::read_dta(paste0(datapath, "out/ecz-BMI-all.dta"))
var_eth <- haven::read_dta(paste0(datapath, "out/ecz-ethnicity.dta"))
var_alc <- haven::read_dta(paste0(datapath, "out/ecz-harmfulalcohol.dta"))
var_smok <- haven::read_dta(paste0(datapath, "out/ecz-smoke-all.dta"))

# get rid of indexdate because they get in the way of merging later
var_bmi <- var_bmi %>% select(-indexdate)
var_smok <- var_smok %>% select(-indexdate)

# merge into one big ol’ dataset ------------------------------------------
df_tab1 <- cohort %>%
	#sample_n(size = 1000) %>%
	left_join(var_agesexgp, by = "patid") %>%
	left_join(var_depressionAll, by = "patid") %>%
	mutate(dep = ifelse(is.na(eventdate), 0, 1)) %>%
	rename(depdate = eventdate) %>%
	left_join(var_anxietyAll, by = "patid") %>%
	mutate(anx = ifelse(is.na(eventdate), 0, 1)) %>%
	rename(anxdate = eventdate) %>%
	left_join(var_bmi, by="patid") %>%
	left_join(var_eth, by="patid") %>%
	left_join(var_alc, by="patid") %>%
	mutate(alc = ifelse(is.na(eventdate), 0, 1)) %>%
	rename(alcdate = eventdate) %>%
	left_join(var_smok, by="patid") %>%
	rename(smokdate = eventdate) 

glimpse(df_tab1)

## get rid of big read in data to save memory
rm(list=c(
	"patinfo",
	"var_depressionDef",
	"var_anxietyDef",
	"var_depressionAll",
	"var_anxietyAll",
	"var_agesexgp",
	"var_bmi",
	"var_eth",
	"var_alc",
	"var_smok"))


# bit of formatting 
df_tab1 <- df_tab1 %>%
	mutate(dob = as.Date(paste(realyob, "01", "01", sep = "-")),
				 age = as.numeric((indexdate - dob)/365.25)) %>%
	mutate_at(c("exposed", "gender", "dep", "anx", "eth5", "eth16", "alc", "smokstatus"), ~as_factor(.))


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
summary(df_tab1)

## skimr
#skim(df_tab1)
group_by(df_tab1, exposed) %>% skim()

## summarytools
#summarytools::descr(df_tab1)
#summarytools::dfSummary(df_tab1)

tab1 <- df_tab1 %>% 
	select(patid, exposed, futime, gender, age, age_cat, eth5, bmi, bmi_cat,
				 dep, anx, alc, smokstatus) %>%
	mutate_at("exposed", ~ifelse(. == 0, "Matched controls", "With eczema"))

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
													 smokstatus = "Smoking status")
	) %>%
	add_overall() %>%
	bold_labels() %>%
	modify_footnote(
		all_stat_cols() ~ "Median (IQR) or Frequency (%)"
	) 

table1

table1 %>%
	as_flex_table() %>%
	flextable::save_as_docx(path = here::here("out/tab1.docx"))



# run an easy model  ------------------------------------------------------
cohort
analysis <- tab1 %>%
	select(-exposed) %>%
	left_join(cohort, by = c("patid")) %>%
	select(-bign) %>%
	mutate(time = as.numeric(enddate-indexdate))
haven::write_dta(data = analysis, path = paste0(datapath,"analysis/eczema_analysis.dta"))

#Merge and create time dependent covariates and events
cohort_obs <- cohort %>% 
	group_by(patid) %>% 
	mutate(obs = 1:n()) 
## get those that are not duplicated
cohort_nodups <- cohort_obs %>%
	filter(obs==1)
## get duplciates and rename patid for exposed period 
cohort_dups <- cohort_obs %>%
	filter(max(obs)>1) %>%
	ungroup() %>%
	filter(exposed==0) %>%
	mutate(newpatid = as.numeric(paste0("9999", patid))) %>%
		select(-patid) %>%
		select(patid = newpatid, everything())
	
cohort_edited <- cohort_nodups %>%
	bind_rows(cohort_dups)
nrow(cohort)
nrow(cohort_edited)
sum(duplicated(cohort$patid))
sum(duplicated(cohort_edited$patid))

cohort_sample <- cohort_edited[1:150,]
anx_tmerged0 <- tmerge(cohort_edited, cohort_edited, id=patid, tstart = indexdate, tstop = enddate) #This first step sets the range of follow up
anx_tmerged1 <- tmerge(tmerged0, var_depressionAll, id=patid, dep=tdc(eventdate))
anx_tmerged2 <- tmerge(tmerged1, var_anxietyAll, id=patid, anx=event(eventdate))
# 
# tmerged2 <- tmerge(tmerged1, prescriptions, id=patid, gc=tdc(start), high=tdc(start, high, 0))
# tmerged3 <- tmerge(tmerged2, presc_exp, id=patid, active=tdc(start, active, "not active"))
# tmerged4 <- tmerge(tmerged3, asthma, id=patid, asthma=tdc(eventdate))
# tmerged5 <- tmerge(tmerged4, alc, id=patid, alc=tdc(eventdate))
# tmerged6 <- tmerge(tmerged5, severity, id=patid, severity=tdc(date, modsevere, 0))
# tmerged7 <- tmerge(tmerged6, allfract[!is.na(allfract$include_either),], id=patid, fract_any=event(eventdate, include_either))
# tmerged8 <- tmerge(tmerged7, allfract, id=patid,
# 									 fract_humerus=event(eventdate, prox_hum),
# 									 fract_pelvis=event(eventdate, Pelvis),
# 									 fract_hip=event(eventdate, Hip),
# 									 fract_spine=event(eventdate, spine),
# 									 fract_wrist=event(eventdate, wrist))


tmerged2$t <- (tmerged2$tstop-tmerged2$tstart)/365.25

KM <- survfit(Surv(t, anx) ~ exposed, data = tmerged2)


ggsurvplot(KM, data = tmerged2,conf.int = T,ylim=c(0,1),censor=F,legend.title="Eczema~anx",legend.labs = c("unexposed","exposed"))
dev.copy(pdf, file = here::here("out/km_anxiety.pdf"), 6,6)
	dev.off()
KMunicate(fit = KM, time_scale = KM$time)

dep_tmerged1 <- tmerge(tmerged0, var_anxietyAll, id=patid, anx=tdc(eventdate))
dep_tmerged2 <- tmerge(tmerged1, var_depressionAll, id=patid, dep=event(eventdate))
dep_tmerged2$t <- (tmerged2$tstop-tmerged2$tstart)/365.25

KM2 <- survfit(Surv(t, dep) ~ exposed, data = dep_tmerged2)

pdf(file = here::here("out/km_depression.pdf"), 6,6)
	ggsurvplot(KM2, data = tmerged2,conf.int = T,ylim=c(0,1),censor=F,legend.title="Eczema~dep",legend.labs = c("unexposed","exposed"))
dev.off()



anx_km <- survfit(Surv(time,anx)~as.factor(exposed), data=analysis)
plot(anx_km,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,0.8,c("None","Exposed"),col=c("black","grey"),lty=1,cex=0.5)

plot(anx_km,mark.time=F,col=c("blue","red"),lty=2,add=F,xlab="Time",ylab="Estimated survivor function")
legend(8,0.8,c("None","Exposed"),col=c("blue","red"),lty=1,cex=0.5)

dev.copy(pdf, file = here::here("out/km_anxiety.pdf"), 6,6)
	dev.off()

with(analysis, table(exposed, anx))

dep_km <- survfit(Surv(time,dep)~as.factor(exposed), data=analysis)
plot(dep_km,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,1,c("None","Exposed"),col=c("black","grey"),lty=1,cex=0.5)
with(analysis, table(exposed, dep))

dep_cox<-coxph(Surv(time,dep)~as.factor(exposed),data=analysis)
summary(dep_cox)

dep_survfit <- survfit(dep_cox,newdata=data.frame(exposed=c(0,1)))
summary(dep_survfit)

plot(dep_survfit,mark.time=F,col=c("black","grey"),xlab="Time",ylab="Estimated survivor function")
legend(8,1,c("Non","Eczema"),col=c("black","grey"),lty=1,cex=0.5)


# Old gt version ----------------------------------------------------------
# 
# # build a table 1 with `gt` -----------------------------------------------
# # Numeric (age, fuptime, BMI)
# numvars <- df_tab1 %>%
# 	group_by(exposed) %>% 
# 	summarise(n = n(),
# 						sumpyears = sum(futime),
# 						medfup = median(futime), 
# 						iqrfupL = quantile(futime, 0.25),
# 						iqrfupU = quantile(futime, 0.75),
# 						medage = median(age), 
# 						iqrageL = quantile(age, 0.25, na.rm = T),
# 						iqrageU = quantile(age, 0.75, na.rm = T),
# 						medbmi = median(bmi, na.rm = T), 
# 						iqrbmiL = quantile(bmi, 0.25, na.rm = T),
# 						iqrbmiU = quantile(bmi, 0.75, na.rm = T)
# 						) 
# # format the table 
# numvars_fmt <- numvars	%>% pivot_longer(-exposed, names_to = "Category") %>%
# 	pivot_wider(names_from = exposed) %>%
# 	rename(n0 = `0`, n1 = `1`) %>%
# 	mutate(var = c("N", 
# 								 rep("Follow-up",4),
# 								 rep("Age",3),
# 								 rep("BMI",3)
# 	)) 
# 
# 
# numvars_iqr <- numvars_fmt %>%
# 	mutate(iqr = ifelse(grepl("iqr", Category) & grepl("L", Category), "iqrL", "n")) %>%
# 	mutate_at("iqr" , ~ifelse(. == "n" & grepl("iqr", Category) & grepl("U", Category), "iqrU", .)) %>%
# 	filter(grepl("iqr", iqr)) %>%
# 	rename(pc0 = n0, pc1 = n1) %>% ## pc is really iqr %>%
# 	select(-Category) %>%
# 	pivot_wider(id_cols = var, names_from = iqr, values_from = c(pc0, pc1)) 
# 
# numvars_fmt2 <- numvars_fmt %>%
# 	filter(!grepl("iqr", Category)) %>%
# 	left_join(numvars_iqr, by = "var") %>%
# 	mutate_at(c("pc0_iqrL", "pc0_iqrU", "pc1_iqrL", "pc1_iqrU"), 
# 						~ifelse(Category == "sumpyears", NA, .)) %>%
# 	mutate_at("Category",
# 						~ifelse(.=="sumpyears", "Total person-years",
# 										ifelse(.=="medfup", "Duration of follow-up, median (IQR), yrs",
# 													 ifelse(.=="medage", "Age, median (IQR), yrs",
# 													 			 ifelse(.=="medbmi", "BMI, median (IQR), yrs", NA
# 										)))))
# 
# 	
# # Categorical  (sex, ethnicity, agecat, BMIcat, alc, smok, Comorbidities)
# cat_fn <- function(var, name){
# 	temp <- df_tab1 %>%
# 		select(var, exposed) %>%
# 		rename(temp = 1)
# 	out <- split(temp, temp$exposed) %>%
# 		purrr::map(tabyl, temp)
# 	out_0 <- out$`0` %>%
# 		rename(Category = temp) %>%
# 		mutate(var = name) %>%
# 		rename(n0 = n, pc0 = percent)
# 	out_1 <- out$`1` %>%
# 		rename(Category = temp) %>%
# 		mutate(var = name) %>%
# 		rename(n1 =n, pc1 = percent)
# 	full_join(out_0, out_1, by = c("Category", "var"))
# }
# 
# anx_cattab <- cat_fn(var = "anx", name = "Anxiety")
# dep_cattab <- cat_fn(var = "dep", name = "Depression")
# age_cattab <- cat_fn(var = "age_cat", name = "Age")
# gendertab <- cat_fn(var = "gender", name = "Gender")
# eth5tab <- cat_fn(var = "eth5", name = "Ethnicity")
# bmi_cattab <- cat_fn(var = "bmi_cat", name = "BMI")
# alctab <- cat_fn(var = "alc", name = "Harmful alcohol use")
# smokstatustab <- cat_fn(var = "smokstatus", name = "Smoking")
# 
# catvars <- bind_rows(
# 	anx_cattab,
# 	dep_cattab,
# 	age_cattab,
# 	gendertab,
# 	eth5tab,
# 	bmi_cattab,
# 	alctab,
# 	smokstatustab
# 	) %>%
# 	select(-c(valid_percent.x, valid_percent.y))
# 
# 
# # make it pretty ----------------------------------------------------------
# tab1 <- numvars_fmt2 %>% 
# 	bind_rows(catvars) %>%
# 	select(var, Category, everything()) %>%
# 	gt(groupname_col = "var") %>%
# 	tab_header(title = "Table 1: Baseline characteristics of the CPRD eczema cohort",
# 						 subtitle = "Data are numbers (%)") %>%
# 	opt_align_table_header(align = "left") %>%
# 	fmt_number(columns = c(3,4), decimals = 1, drop_trailing_zeros = T) %>%
# 	fmt_percent(columns = c(9:10), rows = !is.na(pc0), 
# 							decimals = 1, pattern = "({x})") %>%
# 	cols_merge(columns = c("n0", "pc0_iqrL", "pc0_iqrU"), pattern = "{1} ({2}-{3})") %>%
# 	cols_merge(columns = c("n1", "pc1_iqrL", "pc1_iqrU"), pattern = "{1} ({2}-{3})") %>%
# 	cols_merge(columns = c("n0", "pc0"), pattern = "{1} {2}") %>%
# 	cols_merge(columns = c("n1", "pc1"), pattern = "{1} {2}") %>%
# 	tab_source_note(
# 		source_note = "Abbreviations: BMI, Body Mass Index; IQR, interquartile range; NA, not applicable."
# 	) %>%
# 	tab_footnote(
# 		footnote = "Test",
# 		locations = cells_body(columns = var, rows = Category == "sumpyears")  
# 	)
# tab1
# 
# tab1A <- catvars %>% 
# 	select(var, Category, everything()) %>%
# 	gt(groupname_col = "var") %>%
# 	tab_header(title = "Table 1: Baseline characteristics of the CPRD eczema cohort",
# 						 subtitle = "Data are numbers (%)") %>%
# 	opt_align_table_header(align = "left") %>%
# 	fmt_number(columns = c(3:4), decimals = 1, drop_trailing_zeros = T) %>%
# 	fmt_percent(columns = c(5:6), decimals = 1, pattern = "({x})") %>%
# 	cols_merge(columns = c("n0", "pc0"), pattern = "{1} {2}") %>%
# 	cols_merge(columns = c("n1", "pc1"), pattern = "{1} {2}")
# 	
# tab1B <- numvars_fmt2 %>% 
# 	select(var, Category, everything()) %>%
# 	gt(groupname_col = "var") %>%
# 	tab_header(title = "Table 1: Baseline characteristics of the CPRD eczema cohort",
# 						 subtitle = "Data are numbers (%)") %>%
# 	opt_align_table_header(align = "left") %>%
# 	fmt_number(columns = c(3:8), decimals = 1, drop_trailing_zeros = T) %>%
# 	cols_merge(columns = c("n0", "pc0_iqrL", "pc0_iqrU"), pattern = "{1} ({2}-{3})") %>%
# 	cols_merge(columns = c("n1", "pc1_iqrL", "pc1_iqrU"), pattern = "{1} ({2}-{3})")
# 	
# testA <- tab1A %>% as_raw_html() 
# testB <- tab1B %>% as_raw_html() 
# 
# data.frame(testA) %>% gt() %>% fmt_markdown(columns = everything())
# data.frame(testB) %>% gt() %>% fmt_markdown(columns = everything())
# tab1 %>% gtsave(filename = here::here("out/tab1.rtf"))
# 	
