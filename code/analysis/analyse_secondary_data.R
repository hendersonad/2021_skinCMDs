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


XX <- c("psoriasis", "eczema")
YY <- c("anxiety", "depression")
exposure <- XX[2]
outcome <- YY[1]
df_model_anx <-
	readRDS(paste0(
		datapath,
		"out/models_data/df_model",
		ABBRVexp,
		"_anxiety.rds"
	))

ABBRVexp <- substr(exposure, 1, 3)

df_model <- df_model_anx

mod5_anx <-
	readRDS(paste0(
		datapath,
		"out/models_data/",
		ABBRVexp,
		"_anxiety_mod5_interaction_age_modeldata.rds"
	))

mod5_tab <- mod5_anx %>%
	tbl_regression(exp = T,
								 include = c(exposed, exposed:agegroup, agegroup),
								 conf.level = 0.99) %>%
	add_glance_source_note(include = n)



mod6_anx <-
	readRDS(paste0(
		datapath,
		"out/models_data/",
		ABBRVexp,
		"_anxiety_mod6_interaction_comorbid_modeldata.rds"
	))

st <- Sys.time()
mod6_tab <- mod6_anx %>%
	tbl_regression(exp = T,
								 include = c(exposed, exposed:comorbid, comorbid),
								 conf.level = 0.99)
end <- Sys.time()
end-st



