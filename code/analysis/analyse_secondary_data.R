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

dir.create(file.path(here("out")), showWarnings = FALSE)
dir.create(file.path(here("out", "analysis")), showWarnings = FALSE)


YY <- c("depression", "anxiety")
XX <- c("psoriasis", "eczema")

st <- Sys.time()
for(exposure in XX){
	ABBRVexp <- substr(exposure, 1 , 3)
	
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
	  .dib(paste0(outcome, "~", exposure))
	  # load models -------------------------------------------------------------
	  mod4 <-
	    readRDS(
	      paste0(
	        datapath,
	        "out/models_data/",
	        ABBRVexp,
	        "_",
	        outcome,
	        "_mod4_severity_modeldata.rds"
	      )
	    )
	  mod5 <-
	    readRDS(
	      paste0(
	        datapath,
	        "out/models_data/",
	        ABBRVexp,
	        "_",
	        outcome,
	        "_mod5_interaction_age_modeldata.rds"
	      )
	    )
	  mod6 <-
	    readRDS(
	      paste0(
	        datapath,
	        "out/models_data/",
	        ABBRVexp,
	        "_",
	        outcome,
	        "_mod6_interaction_comorbid_modeldata.rds"
	      )
	    )
	  mod7 <-
	    readRDS(
	      paste0(
	        datapath,
	        "out/models_data/",
	        ABBRVexp,
	        "_",
	        outcome,
	        "_mod7_interaction_calendar_modeldata.rds"
	      )
	    )
	  
	  out <- substr(outcome, 1, 3)
	  df_model <- get(paste0("df_model_", out))
	  
	  mod4_tab <- mod4 %>%
	    tbl_regression(exp = T,
	                   #include = "severity",
	                   conf.level = 0.99)
	  
	  mod5_tab <- mod5 %>%
	    tbl_regression(
	      exp = T,
	      include = c('exposed', 'exposed:agegroup', 'agegroup'),
	      conf.level = 0.99
	    ) 
	  
	  mod6_tab <- mod6 %>%
	    tbl_regression(
	      exp = T,
	      include = c('exposed', 'exposed:comorbid', 'comorbid'),
	      conf.level = 0.99
	    )
	  
	  mod7_tab <- mod7 %>%
	    tbl_regression(
	      exp = T,
	      include = c('exposed', 'exposed:cal_period', 'cal_period'),
	      conf.level = 0.99
	    )
	  
	  # summarise results -------------------------------------------------------
	  tbls_secondary <- tbl_merge(
	    tbls = list(mod5_tab, mod6_tab, mod7_tab),
	    tab_spanner = c("**Age group**",
	                    "**Comorbidity**",
	                    "**Calendar Period**")
	  )
	  
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
	rm(mod4, mod5, mod6, mod7,
	   mod4_tab, mod5_tab, mod6_tab, mod7_tab)
}

end <- Sys.time()
end-st
