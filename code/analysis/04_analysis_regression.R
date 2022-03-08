library(tidyverse)
library(here)
library(labelled)
library(stringr)
library(janitor)
library(timetk)
library(gtsummary)
library(magrittr)
library(survival)
library(survminer)
library(htmltools)
library(lubridate)

if (Sys.info()["user"] == "lsh1510922") {
	if (Sys.info()["sysname"] == "Darwin") {
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if (Sys.info()["sysname"] == "Windows") {
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")
samplingSmall <- F
export_plots <- F # will be written to F if samplingSmall == T
export_models <- T
export_tables <- F
small_sample_size <- 100

st_time <- Sys.time()
for (exposure in XX) {
	#exposure <- XX[1]
	ABBRVexp <- str_sub(exposure, 1 , 3)
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
	
	.dib(exposure)
	
	
	if (samplingSmall) {
		patids_anx <-  sample(unique(df_anx_split$setid), size = small_sample_size)
		patids_dep <-  sample(unique(df_dep_split$setid), size = small_sample_size)
		
		df_anx_split <- df_anx_split %>% filter(setid %in% patids_anx)
		df_dep_split <- df_dep_split %>% filter(setid %in% patids_dep)
		export_plots <- F
	}
	
	df_anx_split$t <-
		as.numeric(df_anx_split$tstop - df_anx_split$tstart)
	df_dep_split$t <-
		as.numeric(df_dep_split$tstop - df_dep_split$tstart)
	
	df_anx_split$gc90days <-
		factor(df_anx_split$gc90days, levels = c("not active", "active"))
	df_dep_split$gc90days <-
		factor(df_dep_split$gc90days, levels = c("not active", "active"))
	
	df_anx_split$bmi2 <-
		df_anx_split$bmi - mean(df_anx_split$bmi, na.rm = T)
	df_dep_split$bmi2 <-
		df_dep_split$bmi - mean(df_dep_split$bmi, na.rm = T)
	
	# annual period prevalence plot -------------------------------------------
	if (export_plots) {
		temp <- df_anx_split %>%
			select(setid, patid, exposed, dob, indexdate, tstart, tstop, out) %>%
			mutate(eligible_from = dob + tstart,
						 eligible_to = dob + tstop)
		
		#Convert tstart and tstop back to dates
		df_split_cal <- survSplit(
			Surv(
				time = as.numeric(eligible_from),
				time2 = as.numeric(eligible_to),
				event = out
			) ~ .,
			data = temp,
			cut = as.numeric(as.Date("1996-01-02") + c(365.25 * 1:23)),
			episode = "year"
		)
		
		
		df_split_cal <- df_split_cal %>%
			group_by(setid, patid, year) %>%
			slice(1)
		df_split_cal$year2 <- factor(df_split_cal$year,
																 levels = 1:24,
																 labels = as.character(1997:2020))
		df_split_sum <- df_split_cal %>%
			group_by(exposed, year2) %>%
			summarise(at_risk = n(), outcome = sum(out)) %>%
			mutate(prop = outcome / at_risk,
						 exposed = factor(
						 	exposed,
						 	levels = 0:1,
						 	labels = c("Matched controls", paste0("With ", exposure))
						 )) %>%
			ungroup() %>%
			rowwise() %>%
			mutate(tst = list(broom::tidy(
				prop.test(outcome, at_risk, conf.level = 0.95)
			))) %>%
			tidyr::unnest(tst)
		
		p1 <-
			ggplot(
				df_split_sum,
				aes(
					x = year2,
					y = estimate * 100,
					ymin = conf.low * 100,
					ymax = conf.high * 100,
					group = exposed,
					colour = exposed,
					fill = exposed
				)
			) +
			geom_line(lty = 2) +
			geom_ribbon(lty = 0, alpha = 0.2) +
			labs(title = paste0(exposure, " ~ anxiety"),
					 y = "Period prevalence",
					 x = "Year") +
			theme_ali() +
			theme(axis.text.x = element_text(angle = 30))
		
		dodge <- position_dodge(width = 0.9)
		p2 <-
			ggplot(
				df_split_sum,
				aes(
					x = year2,
					y = estimate * 100,
					ymin = conf.low * 100,
					ymax = conf.high * 100,
					group = exposed,
					colour = exposed,
					fill = exposed
				)
			) +
			geom_point(pch = 16, position = dodge) +
			geom_errorbar(lty = 1,
										alpha = 0.8,
										position = dodge) +
			labs(title = paste0(exposure, " ~ anxiety"),
					 y = "Period prevalence",
					 x = "Year") +
			theme_ali() +
			theme(axis.text.x = element_text(angle = 30))
		
		
		temp <- df_dep_split %>%
			select(setid, patid, exposed, dob, indexdate, tstart, tstop, out) %>%
			mutate(eligible_from = dob + tstart,
						 eligible_to = dob + tstop)
		
		#Convert tstart and tstop back to dates
		df_split_cal <- survSplit(
			Surv(
				time = as.numeric(eligible_from),
				time2 = as.numeric(eligible_to),
				event = out
			) ~ .,
			data = temp,
			cut = as.numeric(as.Date("1996-01-02") + c(365.25 * 1:23)),
			episode = "year"
		)
		
		
		df_split_cal <- df_split_cal %>%
			group_by(setid, patid, year) %>%
			slice(1)
		df_split_cal$year2 <- factor(df_split_cal$year,
																 levels = 1:24,
																 labels = as.character(1997:2020))
		df_split_sum <- df_split_cal %>%
			group_by(exposed, year2) %>%
			summarise(at_risk = n(), outcome = sum(out)) %>%
			mutate(prop = outcome / at_risk,
						 exposed = factor(
						 	exposed,
						 	levels = 0:1,
						 	labels = c("Matched controls", paste0("With ", exposure))
						 )) %>%
			ungroup() %>%
			rowwise() %>%
			mutate(tst = list(broom::tidy(
				prop.test(outcome, at_risk, conf.level = 0.95)
			))) %>%
			tidyr::unnest(tst)
		
		p3 <-
			ggplot(
				df_split_sum,
				aes(
					x = year2,
					y = estimate * 100,
					ymin = conf.low * 100,
					ymax = conf.high * 100,
					group = exposed,
					colour = exposed,
					fill = exposed
				)
			) +
			geom_line(lty = 2) +
			geom_ribbon(lty = 0, alpha = 0.2) +
			labs(
				title = paste0(exposure, " ~ depression"),
				y = "Period prevalence",
				x = "Year"
			) +
			theme_ali() +
			theme(axis.text.x = element_text(angle = 30))
		
		
		#dodge <- position_dodge(width=0.9)
		p4 <-
			ggplot(
				df_split_sum,
				aes(
					x = year2,
					y = estimate * 100,
					ymin = conf.low * 100,
					ymax = conf.high * 100,
					group = exposed,
					colour = exposed,
					fill = exposed
				)
			) +
			geom_point(pch = 16, position = dodge) +
			geom_errorbar(lty = 1,
										alpha = 0.8,
										position = dodge) +
			labs(
				title = paste0(exposure, " ~ depression"),
				y = "Period prevalence",
				x = "Year"
			) +
			theme_ali() +
			theme(axis.text.x = element_text(angle = 30))
		
		pAll <- cowplot::plot_grid(p2, p4, ncol = 1)
		
		pdf(here::here(
			"out/analysis",
			paste0("ann-prevalence-anxiety-", ABBRVexp, ".pdf")
		), 10, 6)
		print(p1)
		dev.off()
		pdf(here::here(
			"out/analysis",
			paste0("ann-prevalence-anxiety-", ABBRVexp, "_v2.pdf")
		), 10, 6)
		print(p2)
		dev.off()
		pdf(here::here(
			"out/analysis",
			paste0("ann-prevalence-depression-", ABBRVexp, ".pdf")
		), 10, 6)
		print(p3)
		dev.off()
		pdf(here::here(
			"out/analysis",
			paste0("ann-prevalence-depression-", ABBRVexp, "_v2.pdf")
		), 10, 6)
		print(p4)
		dev.off()
		pdf(here::here(
			"out/analysis",
			paste0("ann-prevalence-", ABBRVexp, ".pdf")
		), 10, 8)
		print(pAll)
		dev.off()
		
		rm(df_split_sum, df_split_cal, temp, p1, p2, p3, p4)
		
		# Run a simple Exponential regression --------------------------------
		#.dib("Exp model")
		#
		#survreg(Surv(t) ~ exposed, dist = "exponential", data = df_anx_split) %>%
		#	tbl_regression(exp = TRUE)
		
		# Run a simple KM analysis ------------------------------------------------
		pdf(file = here::here("out", paste0("km_", ABBRVexp, "_anxiety_v2.pdf")),
				6,
				6,
				onefile = F)
		km <- ggsurvplot(
			fit = survminer::surv_fit(Surv(t / 365.25, out) ~ exposed, data = df_anx_split),
			xlab = "Years",
			ylab = "Overall survival probability",
			conf.int = T,
			censor = F,
			ylim = c(0.9, 1),
			legend.title = paste0(exposure, "~anx"),
			legend.labs = c("unexposed", "exposed")
		)
		print(km)
		dev.off()
		### improvements: add in a little tile with c(0,1) ylims to show true survival
		
		pdf(file = here::here("out", paste0(
			"km_", ABBRVexp, "_depression_v2.pdf"
		)),
		6,
		6,
		onefile = F)
		km <- ggsurvplot(
			fit = survminer::surv_fit(Surv(t / 365.25, out) ~ exposed, data = df_dep_split),
			xlab = "Years",
			ylab = "Overall survival probability",
			conf.int = T,
			censor = F,
			ylim = c(0.9, 1),
			legend.title = paste0(exposure, "~anx"),
			legend.labs = c("unexposed", "exposed")
		)
		print(km)
		dev.off()
		rm(km)
	}
	
	for (outcome in YY) {
		#outcome = YY[2]
		.dib(outcome)
		
		if (outcome == "anxiety") {
			df_model <- df_anx_split
		} else if (outcome == "depression") {
			df_model <- df_dep_split
		}

		# Bit of variable formatting for output in regression tables --------------
		df_model <- df_model %>% 
			mutate(exposed = case_when(
				exposed == 0 ~ "Unexposed",
				exposed == 1 ~ stringr::str_to_title(paste0(exposure))
			))
			#rename(!!paste0(exposure) := exposed) %>% 
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
		
		saveRDS(df_model, file = paste0(datapath, "out/models_data/df_model", ABBRVexp, "_", outcome, ".rds"))
		
		# Run a simple Cox regression ----------------------------------------
		.dib("Running model 1 (crude)")
		mod1 <-
			coxph(Surv(t, out) ~ exposed + strata(setid), data = df_model) #%>%
		#mod1_tab <- mod1 %>% 
		#	tbl_regression(exp = TRUE) %>%
		#	add_n() %>% 
		#	modify_footnote(update = estimate ~ "Adjusted for matched set (age, sex, GP)")
		#mod1
		
		# Run a confounder Cox regression ----------------------------------------
		.dib("Running model 2 (confounder)")
		mod2 <-
			coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
						data = df_model)
		# mod2_exp <- mod2 %>% tbl_regression(exp = T,
		# 																		include = "exposed",
		# 																		conf.level = 0.99) %>% 
		# 	add_n() %>% 
		# 	modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
		# mod2_full <- mod2 %>% tbl_regression(exp = T, conf.level = 0.99) %>% 
		# 	modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
		#mod2
		
		# Run a mediator Cox regression -------------------------------------------
		.dib("Running model 3 (mediator)")
		if (ABBRVexp == "ecz") {
			mod3 <-
				coxph(
					Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + gc90days + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod3 <-
				coxph(
					Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
					data = df_model
				) 
		}
		#mod3
		# mod3_exp <- mod3 %>% tbl_regression(exp = T,
		# 																		include = "exposed",
		# 																		conf.level = 0.99) %>% 
		# 	add_n() 
		# mod3_full <- mod3 %>% tbl_regression(exp = T, conf.level = 0.99) %>% 
		# 	modify_footnote(update = estimate ~ "Additionally adjusted for calendar period and comorbidities")
		# if(exposure == "eczema"){
		# 	mod3_exp <- mod3_exp %>%
		# 		modify_footnote(update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status, sleep problems, oral glucocorticoid use")
		# 	mod3_full <- mod3_full %>% 
		# 		modify_footnote(update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status, sleep problems, oral glucocorticoid use")
		# }else if (exposure == "psoriasis"){
		# 	mod3_exp <- mod3_exp %>%
		# 		modify_footnote(update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status")
		# 	mod3_full <- mod3_full %>% 
		# 		modify_footnote(update = estimate ~ "Additionally adjusted for BMI, alcohol misuse, smoking status")
		# }
		
		# summarise results -------------------------------------------------------
		# calculate the total N by exposed
		# tbl_n <- df_model %>% select(exposed) %>%
		# 	tbl_summary(statistic = everything() ~ "{n}") %>%
		# 	modify_header(stat_0 ~ "**N**") %>%
		# 	modify_footnote(everything() ~ NA)
		# 
		# # calculate the event N by exposed
		# tbl_nevent <- df_model %>% filter(out == 1) %>% select(exposed) %>%
		# 	tbl_summary(statistic = everything() ~ "{n}") %>%
		# 	modify_header(stat_0 ~ "**Event N**")  %>%
		# 	modify_footnote(everything() ~ NA)
		# 
		# # calculate the pyars N by exposed
		# tbl_pyars <- df_model %>% select(exposed, t) %>%
		# 	group_by(exposed) %>%
		# 	summarise(pyars = round(sum(t))) %>%
		# 	uncount(pyars) %>%
		# 	tbl_summary(statistic = everything() ~ "{n}") %>% 
		# 	modify_header(stat_0 ~ "**P-years**")  %>%
		# 	modify_footnote(everything() ~ NA)
		# 
		# merge tbls together
		# tbls <-
		# 	tbl_merge(
		# 		list(tbl_n, tbl_nevent, tbl_pyars, mod1_tab, mod2_exp, mod3_exp),
		# 		tab_spanner = c(NA, NA, NA, "**Crude**", "**Confounder**", "**Mediator**")
		# 	) 
		# 
		# tbls_full <- tbl_merge(
		# 	tbls = list(tbl_n, tbl_nevent, tbl_pyars, mod1_tab, mod2_full, mod3_full),
		# 	tab_spanner = c(NA, NA, NA, "**Crude**", "**Confounder**", "**Mediator**")
		# )
		
		if (export_tables) {
			tbls %>%
				gtsummary::as_gt() %>%
				gt::gtsave(
					filename =  paste0("tabls1A_", ABBRVexp, "_", outcome, ".html"),
					path = here::here("out/analysis")
				)
			tbls_full %>%
				gtsummary::as_gt() %>%
				gt::gtsave(
					filename =  paste0("tabls1B_", ABBRVexp, "_", outcome, "_full.html"),
					path = here::here("out/analysis")
				)
		}
		
		if (export_models) {
			saveRDS(
				mod1,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod1_modeldata.rds"
				)
			)
			saveRDS(
				mod2,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod2_modeldata.rds"
				)
			)
			saveRDS(
				mod3,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod3_modeldata.rds"
				)
			)
		}
		
		rm(mod1, mod2, mod3)

		# Severity ----------------------------------------------------------------
		.dib("Running model 4 (severity)")
		if (ABBRVexp == "ecz") {
			mod4 <-
				coxph(
					Surv(t, out) ~ severity + carstairs + cal_period + comorbid + cci + bmi2 + sleep + alc + smokstatus + gc90days + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod4 <-
				coxph(
					Surv(t, out) ~ severity + carstairs + cal_period + comorbid + cci + bmi2 + alc + smokstatus + strata(setid),
					data = df_model
				) 
		}
		# mod4_tab <- mod4 %>%
		# 	tbl_regression(exp = T,
		# 								 include = "severity",
		# 								 conf.level = 0.99) %>%
		# 	add_n() %>%
		# 	add_nevent()
		
		# Interaction -------------------------------------------------------------
		.dib("Running model 5 (interaction age)")
		if (ABBRVexp == "ecz") {
			mod5 <-
				coxph(
					Surv(t, out) ~ exposed + agegroup + exposed * agegroup + carstairs + cal_period + comorbid + cci + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod5 <-
				coxph(
					Surv(t, out) ~ exposed + agegroup + exposed * agegroup + carstairs + cal_period + comorbid + cci + strata(setid),
					data = df_model
				) 
		}
		# mod5_tab <- mod5 %>%
		# 	tbl_regression(exp = T,
		# 								 include = c(exposed, exposed:agegroup, agegroup),
		# 								 conf.level = 0.99) %>%
		# 	add_glance_source_note(include = n) 
		
		.dib("Running model 6 (interaction comorbid)")
		if (ABBRVexp == "ecz") {
			mod6 <-
				coxph(
					Surv(t, out) ~ exposed + comorbid + exposed * comorbid + carstairs + cal_period + comorbid + cci + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod6 <-
				coxph(
					Surv(t, out) ~ exposed + comorbid + exposed * comorbid + carstairs + cal_period + comorbid + cci + strata(setid),
					data = df_model
				) 
		}
		# mod6_tab <- mod6 %>%
		# 	tbl_regression(exp = T,
		# 								 include = c(exposed, exposed:comorbid, comorbid),
		# 								 conf.level = 0.99) 
		
		.dib("Running model 7 (interaction calendar)")
		if (ABBRVexp == "ecz") {
			mod7 <-
				coxph(
					Surv(t, out) ~ exposed + cal_period + exposed * cal_period + carstairs + comorbid + cci + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod7 <-
				coxph(
					Surv(t, out) ~ exposed + cal_period + exposed * cal_period + carstairs + comorbid + cci + strata(setid),
					data = df_model
				)
		}
		# mod7_tab <- mod7 %>%
		# 	tbl_regression(exp = T,
		# 								 include = c(exposed, exposed:cal_period, cal_period),
		# 								 conf.level = 0.99) 
		
		# summarise results -------------------------------------------------------
		# tbls_secondary <- tbl_merge(
		# 	tbls = list(mod5_tab, mod6_tab, mod7_tab),
		# 	tab_spanner = c(
		# 		"**Age group**",
		# 		"**Comorbidity**",
		# 		"**Calendar Period**"
		# 	)
		# )
		# 
		# if (export_tables) {
		# 	mod4_tab %>%
		# 		gtsummary::as_gt() %>%
		# 		gt::gtsave(
		# 			filename =  paste0("tabls2_", ABBRVexp, "_", outcome, "_severity.html"),
		# 			path = here::here("out/analysis")
		# 		)
		# 	tbls_secondary %>%
		# 		gtsummary::as_gt() %>%
		# 		gt::gtsave(
		# 			filename =  paste0("tabls3_", ABBRVexp, "_", outcome, "_interaction.html"),
		# 			path = here::here("out/analysis")
		# 		)
		# }
		if (export_models) {
			saveRDS(
				mod4,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod4_severity_modeldata.rds"
				)
			)
			saveRDS(
				mod5,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod5_interaction_age_modeldata.rds"
				)
			)
			saveRDS(
				mod6,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod6_interaction_comorbid_modeldata.rds"
				)
			)
			saveRDS(
				mod7,
				file = paste0(
					datapath,
					"out/models_data/",
					ABBRVexp,
					"_",
					outcome,
					"_mod7_interaction_calendar_modeldata.rds"
				)
			)
		}
		rm(mod4, mod5, mod6, mod7)
	}
	# merge multiple tables ---------------------------------------------------
	## but they need different names
	
	if (sum(grepl(pattern = "df_model", x = ls())) == 1) {
		rm(df_model)
	}
	if (sum(grepl(pattern = "df_dep_split", x = ls())) == 1) {
		rm(df_dep_split)
	}
	if (sum(grepl(pattern = "df_anx_split", x = ls())) == 1) {
		rm(df_anx_split)
	}
}
end_tim <- Sys.time()

end_tim - st_time
