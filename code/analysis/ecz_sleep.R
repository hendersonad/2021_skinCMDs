# install.packages("rtools")
# install.packages("gt")
# install.packages("gtsummary")
# install.packages("janitor")
# install.packages("here")

library(tidyverse)
library(here)
library(flextable)
library(gt)
library(gtsummary)
library(janitor)
library(timetk)
library(skimr)
library(glue)


if(Sys.info()["user"]=="lsh1510922"){
	if(Sys.info()["sysname"]=="Darwin"){
		datapath <- "/Volumes/EHR Group/GPRD_GOLD/Ali/2021_skinepiextract/"
		datapath <- "/Users/lsh1510922/Documents/Postdoc/2021_extract/"
	}
	if(Sys.info()["sysname"]=="Windows"){
		datapath <- "Z:/GPRD_GOLD/Ali/2021_skinepiextract/"
	}
}

exposure <- "eczema"
ABBRVexp <- str_sub(exposure,1 ,3)

.dib(exposure)

df_anx_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-anxiety_split.rds"))
df_dep_split <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_split.rds"))

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

for (outcome in c("depression", "anxiety")) {
	if(outcome == "depression"){
		df_in = df_dep_split
	}
	if(outcome == "anxiety"){
		df_in = df_anx_split
	}
	df_out <- df_in %>% 
		arrange(setid, patid, tstart) %>% 
		group_by(setid, patid) %>% 
		slice(1)
	
	
	tab_sleep <- twoXtwo(df = df_out, "severity", out = "sleep")
	tab_sleep_all <-
		twoXtwo(df = df_out, "severity", out = "sleep_all")
	
	tab_sleep_out <- bind_rows(tab_sleep, tab_sleep_all)
	tab_sleep_out %>%
		select(-exposure,-Miss) %>%
		gt::gt() %>%
		tab_row_group(label = "Definite sleep",
									rows = 1:5) %>%
		tab_row_group(label = "Probable sleep (incl. benzo)",
									rows = 6:10) %>%
		gt::tab_header(title = "Sleep problems by eczema severity") %>%
		gt::fmt_number(columns = c(3, 5), decimals = 1) %>%
		gt::fmt_number(columns = c(2, 4), decimals = 0) %>%
		gt::data_color(
			columns = c(Yes_pc),
			colors = scales::col_numeric(
				palette = paletteer::paletteer_c(palette = "harrypotter::gryffindor",
																				 n = 100) %>% as.character(),
				domain = c(
					min(tab_sleep_out$Yes_pc, na.rm = T),
					max(tab_sleep_out$Yes_pc, na.rm = T)
				)
			)
		) %>%
		cols_label(val = "Severity",
							 No_pc = "%",
							 Yes_pc = "%") %>%
		gt::gtsave(
			filename =  paste0("eczema_sleep_", outcome,".html"),
			path = here::here("out/analysis")
		)
}