source(here::here("code/packages.R"))
source(here::here("code/file_paths.R"))

dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data"), showWarnings = FALSE)
dir.create(file.path(datapath, "out/models_data_impute"), showWarnings = FALSE)

XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")

export_plots <- TRUE
export_models <- TRUE

st_time <- Sys.time()
for (exposure in XX) {
	ABBRVexp <- str_sub(exposure, 1 , 3)
	if (exposure == "eczema") {
		df_anx_model <-
			readRDS(paste0(datapath, "out/df_modelecz_anxiety.rds"))
		df_dep_model <-
			readRDS(paste0(datapath, "out/df_modelecz_depression.rds"))
	} else if (exposure == "psoriasis") {
		df_anx_model <-
			readRDS(paste0(datapath, "out/df_modelpso_anxiety.rds"))
		df_dep_model <-
			readRDS(paste0(datapath, "out/df_modelpso_depression.rds"))
	}
	
	
	.dib(exposure)
	for (outcome in YY) {
		.dib(outcome)
		if (outcome == "anxiety") {
			df_model <- df_anx_model
		} else if (outcome == "depression") {
			df_model <- df_dep_model
		}
	  ## Edit gender category because it has redundant groups that have no data
	  df_model <- df_model %>% filter(gender %in% c("Male", "Female"))
	  levels(df_model$gender) <- c(NA, "Male", "Female", NA, NA)
	  df_model$gender <- factor(df_model$gender, levels = c("Female", "Male")) ## set Female as baseline
	  
		# Run a simple Cox regression ----------------------------------------
		.dib("Running model 1 (crude)")
		mod1 <-
			coxph(Surv(t, out) ~ exposed + strata(setid), data = df_model) #%>%
		
		# Run a confounder Cox regression ----------------------------------------
		.dib("Running model 2 (confounder)")
		mod2 <-
			coxph(Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + strata(setid),
						data = df_model)
		
		# Run a mediator Cox regression -------------------------------------------
		.dib("Running model 3 (mediator)")
		if (ABBRVexp == "ecz") {
			mod3 <-
				coxph(
					Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + sleep + gc90days + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod3 <-
				coxph(
					Surv(t, out) ~ exposed + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
					data = df_model
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
		
		# Severity ----------------------------------------------------------------
		.dib("Running model 4 (severity)")
		if (ABBRVexp == "ecz") {
			mod4 <-
				coxph(
					Surv(t, out) ~ severity + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + sleep + gc90days + strata(setid),
					data = df_model
				) 
		} else if (ABBRVexp == "pso") {
			mod4 <-
				coxph(
					Surv(t, out) ~ severity + carstairs + cal_period + comorbid + cci + bmi_cat + alc + smokstatus + strata(setid),
					data = df_model
				) 
		}
		
		# Interaction -------------------------------------------------------------
		.dib("Running model 5 (interaction age)")
		mod5 <-
			coxph(
				Surv(t, out) ~ exposed + agegroup + exposed * agegroup + carstairs + cal_period + comorbid + cci + strata(setid),
				data = df_model
			) 
		.dib("Running model 6 (interaction sex)")
		mod6 <-
			coxph(
				Surv(t, out) ~ exposed + gender + exposed * gender + carstairs + cal_period + comorbid + cci + strata(setid),
				data = df_model
			) 
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
					"_mod6_interaction_sex_modeldata.rds"
				)
			)
		}
	}
}
end_tim <- Sys.time()
end_tim - st_time
