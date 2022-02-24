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

cohort <- haven::read_dta(paste0(datapath, "out/getmatchedcohort-psoriasis-main-mhealth.dta"))
patinfo <- haven::read_dta(paste0(datapath, "in/Patient_extract_pso_extract3_1.dta"))
var_depressionDef <- haven::read_dta(paste0(datapath, "out/variables-pso-depression-definite.dta"))
var_anxietyDef <- haven::read_dta(paste0(datapath, "out/variables-pso-anxiety-definite.dta"))
var_depressionAll <- haven::read_dta(paste0(datapath, "out/variables-pso-depression-all.dta"))
var_anxietyAll <- haven::read_dta(paste0(datapath, "out/variables-pso-anxiety-all.dta"))
var_depressionCensor <- haven::read_dta(paste0(datapath, "out/variables-pso-depression-censor.dta"))
var_anxietyCensor <- haven::read_dta(paste0(datapath, "out/variables-pso-anxiety-censor.dta"))
var_agesexgp <- haven::read_dta(paste0(datapath, "out/variables-pso-age-sex-gp", ".dta"))
var_bmi <- haven::read_dta(paste0(datapath, "out/variables-pso-BMI-all.dta"))
var_eth <- haven::read_dta(paste0(datapath, "out/variables-pso-ethnicity.dta"))
var_alc <- haven::read_dta(paste0(datapath, "out/variables-pso-harmfulalcohol.dta"))
var_smok <- haven::read_dta(paste0(datapath, "out/variables-pso-smoke-all.dta"))
var_sev <- haven::read_dta(paste0(datapath, "out/psoriasis-severity.dta"))

# get rid of indexdate because they get in the way of merging later
var_bmi <- var_bmi %>% select(-indexdate)
var_smok <- var_smok %>% select(-indexdate)

# merge into one big ol’ dataset ------------------------------------------
## depression 
df_tab1 <- cohort %>%
	left_join(var_depressionCensor, by = "patid") %>%
	select(-constype, -readcode, -readterm) %>%
	rename(censorDepDate = eventdate) %>%
	mutate_at("censorDepDate", ~ifelse(is.na(.), as.Date("3000-01-01"), .)) %>%
	left_join(var_depressionAll, by = "patid") %>%
	select(-constype, -readcode, -readterm) %>%
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
	select(-constype, -readcode, -readterm) %>%
	rename(censorAnxDate = eventdate) %>%
	mutate_at("censorAnxDate", ~ifelse(is.na(.), as.Date("3000-01-01"), .)) %>%
	left_join(var_anxietyAll, by = "patid") %>%
	select(-constype, -readcode, -readterm) %>%
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
	left_join(var_bmi, by="patid") %>%
	left_join(var_eth, by="patid") %>%
	left_join(var_alc, by="patid") %>%
	select(-ingredient, -prodcode, -bnftext, -src, -medcode, -readterm) %>%
	mutate(alc = ifelse(is.na(eventdate), 0, 1)) %>%
	rename(alcdate = eventdate) %>%
	left_join(var_smok, by="patid") %>%
	rename(smokdate = eventdate) %>%
	left_join(var_sev, by = "patid") %>%
	mutate_at("src_severetreat" , ~ifelse(date_severe < indexdate, 0, .)) %>%
	mutate_at("src_severetreat" , ~ifelse(date_severe > enddate, 0, .)) %>%
	mutate_at("src_severetreat" , ~replace_na(., 0)) 
	#mutate(severe = ifelse(is.na(date_severe), 0, 1))

glimpse(df_tab1)

write_parquet(df_tab1, sink = paste0(datapath, "out/df_tab1_psoriasis.gz.parquet"))
#df_tab1 <- read_parquet(file = paste0(datapath, "out/df_tab1_psoriasis.gz.parquet"))

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
	mutate_at(c("exposed", "gender", "dep", "anx", "eth5", "eth16", "alc", "smokstatus", "src_severetreat"), ~as_factor(.))

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
summary(df_tab1)

## skimr
#skim(df_tab1)
group_by(df_tab1, exposed) %>% skim()

## summarytools
#summarytools::descr(df_tab1)
#summarytools::dfSummary(df_tab1)

tab1 <- df_tab1 %>% 
	select(patid, exposed, futime, gender, age, age_cat, eth5, bmi, bmi_cat,
				 dep, anx, alc, smokstatus, src_severetreat) %>%
	mutate_at("exposed", ~ifelse(. == 0, "Matched controls", "With psoriasis"))

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
				 anx = factor(anx, levels = 0:1, labels = c("No", "Yes")),
				 severe = factor(src_severetreat, levels = 0:3, labels = c("None", "Systemic", "Phototherapy", "Biologic")))

## investigating the unexposed people with severe psoriasis Rx codes
table(tab1$exposed, tab1$src_severetreat)

unexposed_severe <- tab1 %>%
	filter(severe !="None" & exposed == "Matched controls") %>%
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
	count(obs) ## so almost all are poeple who contribute some unexposed time as a control, then become 
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
	select(-patid, -src_severetreat) %>%
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
	flextable::save_as_docx(path = here::here("out/tab1_pso.docx"))

# run an easy model  ------------------------------------------------------
cohort
analysis <- tab1 %>%
	select(-exposed) %>%
	left_join(cohort, by = c("patid")) %>%
	select(-bign) %>%
	mutate(time = as.numeric(enddate-indexdate))
haven::write_dta(data = analysis, path = paste0(datapath,"analysis/psoriasis_analysis.dta"))

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

## Store the edited cohort sample 
saveRDS(cohort_edited, file = paste0(datapath, "out/pso_cohort_edited.rds"))
if(sum(grepl(pattern = "cohort_edited", x = ls()))==0){
	cohort_edited <- read_rds(paste0(datapath, "out/pso_cohort_edited.rds"))
}
cohort_sample <- cohort_edited[1:150,]

## Use tmerge to build survival dataset
tmerged0 <- tmerge(cohort_edited, cohort_edited, id=patid, tstart = indexdate, tstop = enddate) #This first step sets the range of follow up
anx_tmerged1 <- tmerge(tmerged0, var_depressionAll, id=patid, dep=tdc(eventdate))
anx_tmerged2 <- tmerge(anx_tmerged1, var_anxietyAll, id=patid, anx=event(eventdate))

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


anx_tmerged2$t <- (anx_tmerged2$tstop-anx_tmerged2$tstart)/365.25

KM <- survfit(Surv(t, anx) ~ exposed, data = anx_tmerged2)


pdf(file = here::here("out/km_pso_anxiety.pdf"), 6,6)
	ggsurvplot(KM, data = anx_tmerged2,conf.int = T,ylim=c(0,1),censor=F,legend.title="psoriasis~anx",legend.labs = c("unexposed","exposed"))
dev.off()
#KMunicate(fit = KM, time_scale = KM$time)

dep_tmerged1 <- tmerge(tmerged0, var_anxietyAll, id=patid, anx=tdc(eventdate))
dep_tmerged2 <- tmerge(anx_tmerged2, var_depressionAll, id=patid, dep=event(eventdate))
dep_tmerged2$t <- (dep_tmerged2$tstop-dep_tmerged2$tstart)/365.25

KM2 <- survfit(Surv(t, dep) ~ exposed, data = dep_tmerged2)

pdf(file = here::here("out/km_pso_depression.pdf"), 6,6)
	ggsurvplot(KM2, data = dep_tmerged2,conf.int = T,ylim=c(0,1),censor=F,legend.title="psoriasis~dep",legend.labs = c("unexposed","exposed"))
dev.off()

# Who gets depression?  ---------------------------------------------------
dep_sample <- dep_tmerged2[1:10000, ]
dep_sample <- dep_tmerged2 %>%
	left_join(var_agesexgp, by = "patid") %>%
	left_join(var_bmi, by="patid") %>%
	left_join(var_eth, by="patid") %>%
	left_join(var_alc, by="patid") %>%
	select(-ingredient, -prodcode, -bnftext, -src, -medcode, -readterm) %>%
	mutate(alc = ifelse(is.na(eventdate), 0, 1)) %>%
	rename(alcdate = eventdate) %>%
	left_join(var_smok, by="patid") %>%
	rename(smokdate = eventdate) 

# Set age cats
dep_sample <- dep_sample %>%
	mutate(
		dob = as.Date(paste(realyob,"06","01", sep="-")),
		age = as.numeric((indexdate - dob)/365.25),
		age_cat = as.factor(case_when(
		age < 30 ~ "18 - 29",
		age < 40 ~ "30 - 39",
		age < 50 ~ "40 - 49",
		age < 66 ~ "50 - 65",
		age >= 66 ~ "66+")
	),
	bmi_cat = case_when(
		bmi < 20 ~ "<20",
		bmi < 25 ~ "20-24",
		bmi < 30 ~ "25-29",
		bmi < 66 ~ "30+"))
	
dt_dep = data.table(dep_sample)
dt_dep[, age_x := as.numeric(age_cat)]
ages = levels(dt_dep$age_cat)

dt_dep[, Age2 := as.character(age_cat)]
ages2 = unique(dt_dep[order(age)]$Age2)
dt_dep[, Age2 := factor(Age2, levels = ages2)]
dt_dep[, age2_x := as.numeric(Age2)]

dt_dep[, time := as.numeric(t)]

d_standard = dt_dep[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x)]
d_standard[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
d_standard[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]

# Deaths by "what"
dc = dt_dep[!is.na(exposed)]


dc[, what := factor(ifelse(gender==1, "Male", ifelse(gender==2, "Female", "Other/not-specified")), levels = c("Male", "Female", "Other/not-specified"))]
dc[, what2 := bmi_cat]
dc[, what3 := as.factor(eth5)]
dc[, what4 := as.factor(alc)]

dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, "Pso", "Control"), levels = c("Pso", "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

## WHAT
pdf(file = here::here("out/desc_pso_dep_gender.pdf"), 6,6)
	ggplot() +
		geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
																		 ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
		geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
																				y = rate, yend = rate), size = 0.25, linetype = "22") +
		geom_pointrange(data = dp, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
																	 colour = what, shape = z_exp_label, linetype = z_exp_label), size = 1.4, fatten = 2.4, position = position_dodge(width = 1), stroke = 0.5) +
		scale_x_continuous(breaks = 1:length(ages2), labels = ages2, expand = expansion(0.01)) +
		#scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100), limits = c(0.03, 150), oob = scales::oob_squish) +
		scale_linetype_manual(values = ltys) +
		scale_colour_manual(values = clrs) +
		scale_shape_manual(values = c("Pso" = 5, "Control" = 20)) +
		labs(x =  "Age", 
				 y = "Depression per 100\nyrs of followup" , 
				 colour = "gender", shape = NULL) +
		theme_bw() +
		theme(legend.position = c(0.015, 1.025), legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
					legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
					legend.key.size = unit(0.225, "cm")) + 
		guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
					 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
					 linetype = FALSE)
dev.off()

## WHAT2
dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what2, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, "Pso", "Control"), levels = c("Pso", "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

pdf(file = here::here("out/desc_pso_dep_BMI.pdf"), 6,6)
	ggplot() +
		geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
																		 ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
		geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
																				y = rate, yend = rate), size = 0.25, linetype = "22") +
		geom_pointrange(data = dp, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
																	 colour = what2, shape = z_exp_label, linetype = z_exp_label), size = 1.4, fatten = 2.4, position = position_dodge(width = 1), stroke = 0.5) +
		scale_x_continuous(breaks = 1:length(ages2), labels = ages2, expand = expansion(0.01)) +
		#scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100), limits = c(0.03, 150), oob = scales::oob_squish) +
		scale_linetype_manual(values = ltys) +
		scale_colour_manual(values = clrs) +
		scale_shape_manual(values = c("Pso" = 5, "Control" = 20)) +
		labs(x =  "Age", 
				 y = "Depression per 100\nyrs of followup" , 
				 colour = "BMI", shape = NULL) +
		theme_bw() +
		theme(legend.position = c(0.015, 1.025), legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
					legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
					legend.key.size = unit(0.225, "cm")) + 
		guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
					 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
					 linetype = FALSE)
dev.off()
## WHAT3
dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what3, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, "Pso", "Control"), levels = c("Pso", "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

pdf(file = here::here("out/desc_pso_dep_eth.pdf"), 6,6)
	ggplot() +
		geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
																		 ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
		geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
																				y = rate, yend = rate), size = 0.25, linetype = "22") +
		geom_pointrange(data = dp, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
																	 colour = what3, shape = z_exp_label, linetype = z_exp_label), size = 1.4, fatten = 2.4, position = position_dodge(width = 1), stroke = 0.5) +
		scale_x_continuous(breaks = 1:length(ages2), labels = ages2, expand = expansion(0.01)) +
		#scale_y_log10(breaks = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100), limits = c(0.03, 150), oob = scales::oob_squish) +
		scale_linetype_manual(values = ltys) +
		scale_colour_manual(values = clrs) +
		scale_shape_manual(values = c("Pso" = 5, "Control" = 20)) +
		labs(x =  "Age", 
				 y = "Depression per 100\nyrs of followup" , 
				 colour = "Eth5", shape = NULL) +
		theme_bw() +
		theme(legend.position = c(0.015, 1.025), legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
					legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
					legend.key.size = unit(0.225, "cm")) + 
		guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
					 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
					 linetype = FALSE)
dev.off()