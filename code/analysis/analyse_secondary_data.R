library(tidyverse)
library(here)
library(magrittr)
library(gt)
library(gtsummary)
library(survival)

if (Sys.info()["user"] == "lsh1510922") {
	if (Sys.info()["sysname"] == "Darwin") {
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if (Sys.info()["sysname"] == "Windows") {
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

make_regression_tab <- function(exposure) {
	ABBRVexp <- str_sub(exposure, 1 , 3)
	YY <- c("anxiety", "depression")
	
	# load models -------------------------------------------------------------
	mod4_anx <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_anxiety_mod4_severity_modeldata.rds"
			)
		)
	mod5_anx <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_anxiety_mod5_interaction_age_modeldata.rds"
			)
		)
	mod6_anx <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_anxiety_mod6_interaction_comorbid_modeldata.rds"
			)
		)
	mod7_anx <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_anxiety_mod7_interaction_comorbid_modeldata"
			)
		)
	
	# load models -------------------------------------------------------------
	mod4_dep <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_depression_mod4_severity_modeldata.rds"
			)
		)
	mod5_dep <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_depression_mod5_interaction_age_modeldata.rds"
			)
		)
	mod6_dep <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_depression_mod6_interaction_comorbid_modeldata.rds"
			)
		)
	mod7_dep <-
		readRDS(
			paste0(
				datapath,
				"out/models_data/",
				ABBRVexp,
				"_depression_mod7_interaction_comorbid_modeldata"
			)
		)
	
	# load data ---------------------------------------------------------------
	df_model_anx <-
		readRDS(paste0(
			datapath,
			"out/models_data/df_model",
			ABBRVexp,
			"_anxiety.rds"
		))
	df_model_dep <-
		readRDS(paste0(
			datapath,
			"out/models_data/df_model",
			ABBRVexp,
			"_depression.rds"
		))
	
	for (outcome in YY) {
		out <- substr(outcome, 1,3)
		df_model <- get(paste0("df_model_", out))
		
		mod4 <- get(paste0("mod4_", out))
		mod5 <- get(paste0("mod5_", out))
		mod6 <- get(paste0("mod6_", out))
		mod7 <- get(paste0("mod7_", out))
		
		mod4_tab <- mod4 %>%
			tbl_regression(exp = T,
										 #include = "severity",
										 conf.level = 0.99) %>%
			add_n() %>%
			add_nevent()

		mod5_tab <- mod5 %>%
			tbl_regression(
				exp = T,
				include = c(exposed, exposed:agegroup, agegroup),
				conf.level = 0.99
			) %>%
			add_glance_source_note(include = n)
		
		mod6_tab <- mod6 %>%
			tbl_regression(
				exp = T,
				include = c(exposed, exposed:comorbid, comorbid),
				conf.level = 0.99
			)
		
		mod7_tab <- mod7 %>%
			tbl_regression(
				exp = T,
				include = c(exposed, exposed:cal_period, cal_period),
				conf.level = 0.99
			)
		
		# summarise results -------------------------------------------------------
		tbls_secondary <- tbl_merge(
			tbls = list(mod5_tab, mod6_tab, mod7_tab),
			tab_spanner = c("**Age group**",
											"**Comorbidity**",
											"**Calendar Period**")
		)
		
		if (export_tables) {
			mod4_tab %>%
				gtsummary::as_gt() %>%
				gt::gtsave(
					filename =  paste0("tabls2_", ABBRVexp, "_", outcome, "_severity.html"),
					path = here::here("out/analysis")
				)
			tbls_secondary %>%
				gtsummary::as_gt() %>%
				gt::gtsave(
					filename =  paste0("tabls3_", ABBRVexp, "_", outcome, "_interaction.html"),
					path = here::here("out/analysis")
				)
		}
	}
}

st <- Sys.time()
	make_regression_tab(eczema)
	make_regression_tab(depression)
end <- Sys.time()
end-st
