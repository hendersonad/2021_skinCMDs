# remove.packages("data.table")
# install.packages("data.table", type = "source",
# 								 repos = "https://Rdatatable.gitlab.io/data.table")
library(data.table)
library(tidyverse)
library(here)
library(janitor)
library(timetk)
library(gtsummary)
library(magrittr)
library(survival)
library(timereg)
library(survminer)
library(htmltools)
library(lubridate)
library(mice) # for nelson aalen cum hazard calc
library(broom)
library(stargazer)
library(arm)
library(gt)
library(naniar)

if (Sys.info()["user"] == "lsh1510922") {
	if (Sys.info()["sysname"] == "Darwin") {
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if (Sys.info()["sysname"] == "Windows") {
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")
samplingSmall <- F
export_plots <- T
smpl_small_size <- 1000

#for(exposure in XX) {
exposure <- XX[2]
outcome = YY[2]
ABBRVexp <- str_sub(exposure, 1 , 3)
if (exposure == "eczema") {
	df_anx <-
		readRDS(paste0(datapath, "out/ecz_an_anxiety.rds"))
	df_dep <-
		readRDS(paste0(datapath, "out/ecz_an_depression.rds"))
} else if (exposure == "psoriasis") {
	df_anx <-
		readRDS(paste0(datapath, "out/pso_an_anxiety.rds"))
	df_dep <-
		readRDS(paste0(datapath, "out/pso_an_depression.rds"))
}

.dib(exposure)

if (exposure == "eczema") {
	df_anx_split <-
		readRDS(paste0(datapath, "out/ecz-anxiety_split.rds"))
	df_dep_split <-
		readRDS(paste0(datapath, "out/ecz-depression_split.rds"))
} else if (exposure == "psoriasis") {
	df_anx_split <-
		readRDS(paste0(datapath, "out/pso-anxiety_split.rds"))
	df_dep_split <-
		readRDS(paste0(datapath, "out/pso-depression_split.rds"))
}


if (samplingSmall) {
	set.seed(015235846)
	patids_anx <-
		sample(unique(df_anx_split$setid), size = smpl_small_size)
	patids_dep <-
		sample(unique(df_dep_split$setid), size = smpl_small_size)
	
	df_anx_split <- df_anx_split %>% filter(setid %in% patids_anx)
	df_dep_split <- df_dep_split %>% filter(setid %in% patids_dep)
	export_plots <- F
}

df_anx_split$t <-
	as.numeric(df_anx_split$tstop - df_anx_split$tstart) / 365.25
df_dep_split$t <-
	as.numeric(df_dep_split$tstop - df_dep_split$tstart) / 365.25

df_anx_split$gc90days <-
	factor(df_anx_split$gc90days, levels = c("not active", "active"))
df_dep_split$gc90days <-
	factor(df_dep_split$gc90days, levels = c("not active", "active"))

df_anx_split$bmi2 <-
	df_anx_split$bmi - mean(df_anx_split$bmi, na.rm = T)
df_dep_split$bmi2 <-
	df_dep_split$bmi - mean(df_dep_split$bmi, na.rm = T)



if (outcome == "anxiety") {
	df_model <- df_anx_split
} else if (outcome == "depression") {
	df_model <- df_dep_split
}
if (ABBRVexp == "ecz") {
	df_model$comorbid <- df_model$asthma
} else{
	df_model$comorbid <- df_model$arthritis
}

# collapse the data  ------------------------------------------------------
dt_flat <- setDT(df_model)
# outcome 
dt_flat[, out2 := out[.N], by = list(setid, patid) ]

# sleep
dt_flat[, sleep2 := sleep[1], by = list(setid, patid)]
dt_flat[, sleep2_inc := max(sleep), by = list(setid, patid)]

# alcohol
dt_flat[, alc_numeric := as.numeric(alc)-1, ]
dt_flat[, alc2 := alc_numeric[1], by = list(setid, patid)]
dt_flat[, alc2_inc := max(alc_numeric), by = list(setid, patid)]

# comorbid
dt_flat[, comorbid_numeric := as.numeric(comorbid)-1, ]
dt_flat[, comorbid2 := comorbid_numeric[1], by = list(setid, patid)]
dt_flat[, comorbid2_inc := max(comorbid_numeric), by = list(setid, patid)]

# gc90days
dt_flat[, gc90days2 := gc90days[1], by = list(setid, patid)]
dt_flat[, gc90days2_inc := if(any(gc90days=="active")) "active" else "not active", by = list(setid, patid)]

# smoking
dt_flat[, smok2 := smokstatus[1], by = list(setid, patid)]
dt_flat[, smok2_inc := if(any(smokstatus %in% c("current smoker", "current or ex-smoker"))) "smoker" else "non-smoker", by = list(setid, patid)]

# fup
dt_flat[, fup := as.numeric(enddate-indexdate)/365.25]

select_vars <- c(
		"setid",
		"patid",
		"exposed",
		"out2",
		"fup",
		"indexdate",
		"enddate",
		"dob",
		"sleep2",
		"sleep2_inc",
		"alc2",
		"alc2_inc",
		"comorbid2",
		"comorbid2_inc",
		"gc90days2",
		"gc90days2_inc",
		"smok2",
		"smok2_inc",
		"bmi",
		"bmi2",
		"bmi_cat",
		"carstairs",
		"cci",
		"eth_edited",
		"cal_period",
		"gender",
		"pracid"
	) 

dt_flat2 <- dt_flat[, ..select_vars]
dt_flat2[, sleep_med := ifelse(sleep2_inc == 1 & sleep2 == 0, 1, 0)]
dt_flat2[, alc_med := ifelse(alc2_inc == 1 & alc2 == 0, 1, 0)]
dt_flat2[, comorbid_med := ifelse(comorbid2_inc == 1 & comorbid2 == 0, 1, 0)]
dt_flat2[, gc90days_med := ifelse(gc90days2_inc == "active" & gc90days2 == "not active", 1, 0)]
dt_flat2[, smok3 := as.factor(ifelse(smok2 == "ex-smoker", "ex-smoker", smok2_inc))]
dt_flat2[, exposed := as.factor(exposed)]
dt_flat2[, cal_period := as.factor(cal_period)]
dt_flat2[, age_at_match := as.numeric(indexdate - dob) / 365.25]


## select only one row per patid
dt_flat3 <- dt_flat2[ , .SD[1], by = list(setid, patid)]
setnames(dt_flat3, "out2", "out")

saveRDS(
	dt_flat3,
	file = paste0(
		datapath,
		"out/",
		ABBRVexp,
		"_",
		outcome,
		"_DT.rds"
	)
)
dt_flat3 <- readRDS(
	file = paste0(
		datapath,
		"out/",
		ABBRVexp,
		"_",
		outcome,
		"_DT.rds"
	)
)
tab1 <- table(dt_flat3$exposed, dt_flat3$out)
tab1
prop.table(tab1, 1) %>% signif(digits = 2)

dt_flat3 %>%
	group_by(exposed) %>%
	summarise(mean(fup), sd(fup))
ggplot(dt_flat3, aes(x = fup, group = exposed, fill = exposed)) +
	geom_histogram(bins = 50) +
	theme_ali()



# table 1 -----------------------------------------------------------------
table1 <- dt_flat3 %>% 
	dplyr::select(exposed, out, fup, gender, age_at_match, sleep_med, alc_med, comorbid2, gc90days_med, smok3, bmi, carstairs, cci, eth_edited) %>%
	tbl_summary(by = exposed,
							statistic = list(all_continuous() ~ "{p50} ({p25}-{p75})",
															 all_categorical() ~ "{n} ({p}%)"),
							digits = all_continuous() ~ 1,
							label = list(gender = "Sex",
													 fup = "Follow-up time (years)",
													 cal_period = "Calendar period",
													 age_at_match = "Age",
													 eth_edited = "Ethnicity", 
													 bmi = "BMI",
													 alc_med = "Harmful alcohol use", 
													 comorbid2 = ifelse(exposure=="eczema","Asthma diagnosis","Arthritis diagnosis"),
													 sleep_med = "Sleep problems", 
													 smok3 = "Smoking status", 
													 cci = "Charlson's comorbidity index",
													 carstairs = "Carstairs deprivation quintile",
													 ruc = "Rural/Urban"
							)
	) %>%
	add_overall() %>%
	bold_labels() %>%
	modify_table_styling(align = "right", columns=6:8) %>%
	modify_footnote(
		all_stat_cols() ~ "Median (IQR) or Frequency (%)"
	) 

table1
table1 %>%
	as_gt() %>%
	gt::gtsave(filename = paste0("tab1_",ABBRVexp,"LOGITtest.html"), path = here::here("out"))


# logit with matched vars -------------------------------------------------
df_glm <- dt_flat3 %>% slice(1:1000)
glm1 <-
	glm(out ~ exposed + gender + factor(pracid) + age_at_match,
			family = "binomial",
			data = df_glm) 
summary(glm1)
## Just an obsene number of covariates so ignore pracid for now
df_glm <- dt_flat3

# glm1 <-
# 	glm(out ~ exposed + gender + age_at_match,
# 			family = "binomial",
# 			data = df_glm) 
# summary(glm1)
# summ(glm1, confint = T)
# #tbl_regression(glm1, exp = T)
# 
# glm2 <-
# 	glm(out ~ exposed + gender + age_at_match + carstairs + comorbid2 + bmi2 + cci + cal_period + factor(smokstatus),
# 			family = "binomial",
# 			data = df_glm) 
# summary(glm2)
# # tbl_regression(glm2, exp = T)
# 
# glm3 <-
# 	glm(out ~ exposed + gender + age_at_match + carstairs + comorbid2 + bmi2 + cci + cal_period + factor(smokstatus) + factor(sleep_med),
# 			family = "binomial",
# 			data = df_glm) 
# summary(glm3)
# 
# glm4 <-
# 	glm(out ~ exposed + gender + age_at_match + carstairs + comorbid2 + bmi2 + cci + cal_period + factor(smokstatus) + factor(sleep_med) + factor(alc_med) + factor(gc90days_med),
# 			family = "binomial",
# 			data = df_glm) 
# summary(glm4)
# # 
# # stargazer(
# # 	type = "html",
# # 	out = "out/analysis/glm2.html",
# # 	glm1,
# # 	glm2,
# # 	glm3, 
# # 	glm4,
# # 	column.labels = c("Crude", "Confounder", "Sleep", "All"),
# # 	apply.coef = exp, 
# # 	apply.ci = exp, 
# # 	digits = 2,
# # 	ci = T,
# # 	single.row = TRUE,
# # 	star.cutoffs = c(0.01, 0.001, 0.0001)
# # )
# 
# output_logistic <- function(model) {
# 	CI <- exp(confint.default(model)) %>% as.matrix()
# 	broom::tidy(model) %>%
# 		mutate(or = exp(estimate),
# 					 or_se = signif(or * arm::se.coef(model), 2),
# 					 lci = CI[,1],
# 					 uci = CI[,2],
# 					 p = signif(p.value, 1)) %>%
# 		dplyr::select(term, or, p.value, lci, uci)
# }
# 
# exp_glm1 <- output_logistic(model = glm1) #%>% gt::gt()
# exp_glm2 <- output_logistic(model = glm2) #%>% gt::gt()
# exp_glm3 <- output_logistic(model = glm3) #%>% gt::gt()
# exp_glm4 <- output_logistic(model = glm4) #%>% gt::gt()
# nobs(glm1); nobs(glm2); nobs(glm3); nobs(glm4)
# exp_glm1
# exp_glm2
# exp_glm3
# exp_glm4
# tbl_merge(tbls = list(tbl_regression(glm1, exp =T), tbl_regression(glm2, exp =T)))

pdf(paste0(here("out/analysis/"), ABBRVexp,"_",outcome,"_missing.pdf"), 6,6)
	gg_miss_var(df_glm, show_pct = T)
dev.off()

# clogit ------------------------------------------------------------------
cglm_crude <-
	clogit(out ~ exposed + strata(setid),
				 data = df_glm)
#tbl_regression(cglm, exp = T)
summary(cglm_crude)

cglm <-
	clogit(out ~ exposed + carstairs + comorbid2 + bmi2 + cci + strata(setid),
				 data = df_glm)
#tbl_regression(cglm, exp = T)
summary(cglm)

cglm_med <-
	clogit(out ~ exposed + carstairs + comorbid2 + bmi2 + cci + sleep_med + strata(setid) ,
				 data = df_glm)
#tbl_regression(cglm_med, exp = T)
summary(cglm_med)

cglm_full <-
	clogit(out ~ exposed + carstairs + comorbid2 + bmi2 + cci + smok3 + factor(sleep_med) + factor(alc_med) + factor(gc90days_med) + strata(setid) ,
				 data = df_glm)
#tbl_regression(cglm_med, exp = T)
summary(cglm_full)




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
twoXtwo(df_glm, "exposed", "out")
twoXtwo(df_glm, "sleep_med", "out")
twoXtwo(df_glm, "sleep2", "out")
twoXtwo(df_glm, "sleep2_inc", "out")

twoXtwo(df_glm, "gc90days_med", "out")
twoXtwo(df_glm, "alc_med", "out")
twoXtwo(df_glm, "alc2", "out")
twoXtwo(df_glm, "alc2_inc", "out")

twoXtwo(df_glm, "alc_med", "sleep_med")
twoXtwo(df_glm, "gc90days_med", "sleep_med")


tab_sleep <- twoXtwo(df = df_out, "severity", out = "sleep") 
tab_sleep_all <- twoXtwo(df = df_out, "severity", out = "sleep_all") 
mypal <- paletteer::paletteer_c(palette = "harrypotter::gryffindor", n = 100)
mypal2 <- paletteer::paletteer_c(palette = "harrypotter::slytherin", n = 100)

tab_sleep_out <- bind_rows(tab_sleep, tab_sleep_all)
tab_sleep_out %>% 
	select(-exposure, -Miss) %>% 
	gt::gt() %>% 
	tab_row_group(
		label = "Definite sleep",
		rows = 1:5
	) %>% 
	tab_row_group(
		label = "Probable sleep (incl. benzo)",
		rows = 6:10
	) %>% 
	gt::tab_header(title = "Sleep problems by eczema severity") %>% 
	gt::fmt_number(columns = c(3,5), decimals = 1) %>%
	gt::fmt_number(columns = c(2,4), decimals = 0) %>%
	gt::data_color(
		columns = c(Yes_pc), 
		 colors = scales::col_numeric(
		 	palette = paletteer::paletteer_c(
		 		palette = "harrypotter::gryffindor",
		 		n = 100
		 	) %>% as.character(),
		 	domain = c(min(tab_sleep_out$Yes_pc, na.rm = T), max(tab_sleep_out$Yes_pc, na.rm = T)))
	) %>%  
	cols_label(
		val = "Severity",
		No_pc = "%",
		Yes_pc = "%"
	) %>% 
	gt::gtsave(
		filename =  paste0("eczema_sleep.html"),
		path = here::here("out/analysis")
	)
