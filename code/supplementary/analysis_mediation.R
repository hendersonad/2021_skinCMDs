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

#for(exposure in XX) {
	exposure <- XX[2]
	ABBRVexp <- str_sub(exposure, 1 , 3)
	if (exposure == "eczema") {
		df_anx_split <- readRDS(paste0(datapath, "out/ecz-anxiety_split.rds"))
		df_dep_split <- readRDS(paste0(datapath, "out/ecz-depression_split.rds"))
	} else if (exposure == "psoriasis") {
		df_anx_split <- readRDS(paste0(datapath, "out/pso-anxiety_split.rds"))
		df_dep_split <- readRDS(paste0(datapath, "out/pso-depression_split.rds"))
	}
	
	.dib(exposure)
	
	
	if (samplingSmall) {
		set.seed(015235846)
		patids_anx <-  sample(unique(df_anx_split$setid), size = 1e5)
		patids_dep <-  sample(unique(df_dep_split$setid), size = 1e5)
		
		df_anx_split <- df_anx_split %>% filter(setid %in% patids_anx)
		df_dep_split <- df_dep_split %>% filter(setid %in% patids_dep)
		export_plots <- F
	}
	
	df_anx_split$t <-
		as.numeric(df_anx_split$tstop - df_anx_split$tstart)/365.25
	df_dep_split$t <-
		as.numeric(df_dep_split$tstop - df_dep_split$tstart)/365.25
	
	df_anx_split$gc90days <-
		factor(df_anx_split$gc90days, levels = c("not active", "active"))
	df_dep_split$gc90days <-
		factor(df_dep_split$gc90days, levels = c("not active", "active"))
	
	df_anx_split$bmi2 <-
		df_anx_split$bmi - mean(df_anx_split$bmi, na.rm = T)
	df_dep_split$bmi2 <-
		df_dep_split$bmi - mean(df_dep_split$bmi, na.rm = T)
	
	
	outcome = YY[2]
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
	
	# KM plot to see the general picture--------------------------------------
	km <- survfit(Surv(t, out) ~ exposed, data = df_model)
	
	par(mfrow=c(1,1))
	plot(km, int = F, col = c(2, 4), mark.time = F, xlab = "time", ylab = "Cum Haz KM", fun = "cumhaz", lty = 2)
	
	# nelson-aalen cum hazard -------------------------------------------------
	ch_0 <- nelsonaalen(data = filter(df_model, exposed == 0), 
										timevar = t, 
										statusvar = out)
	ch_1 <- nelsonaalen(data = filter(df_model, exposed == 1), 
											timevar = t, 
											statusvar = out)
	df_model_exp <- df_model %>% filter(exposed == 1) 
	df_model_un <- df_model %>% filter(exposed == 0) 
	points(x = df_model_exp$t, y = ch_1, type = "p", col = 4)
	points(x = df_model_un$t, y = ch_0, col = 2)

	# Nelson-Aalen additive hazard model approach-----------------------------
	df_test <- df_model %>%
		filter(!is.na(carstairs)) %>%
		mutate(id2 = paste0(setid, patid),
					 exposed = factor(exposed))

	## cox based mediation analysis
	mod_confound <- coxph(
			Surv(t, out) ~ exposed + carstairs,
			data = df_test)
	mod_mediate <- coxph(
		Surv(t, out) ~ exposed + carstairs + factor(sleep),
		data = df_test
	)
	summary(mod_confound)		
	summary(mod_mediate)
	
	cox <- tibble(
		var = names(mod_mediate$coefficients),
		beta_con = c(exp(mod_confound$coefficients), NA),
		beta_med = exp(mod_mediate$coefficients)
	) %>%
		mutate(
			change = beta_con - beta_med,
			change_pc = (change/beta_con)*100
		)
	cox ## so lowest deprivation group have 1.3 times greater hazard than group 1. When including sleep this changes to 1.25 suggesting some mediation
	## or: exposed had 1.29 times the hazard but attenuated to 1.26
	
	## model on mediator
	# regression of the mediator (sleep) on the exposure (exposed) adjusted for other confounders 
	mediator_model <- aalen(Surv(t, sleep) ~ const(exposed) + const(carstairs) ,
													data = df_test, 
													id = df_test$patid)
	summary(mediator_model) ## being exposed decreases hazard of sleep problems y 0.5 per 1,000 p-years when adjusted for other stuff
	
	test <- glm(sleep ~ exposed + carstairs , family = "binomial", data = df_test)
	summary(test) ## exposed increases log(odds) of sleep med by 0.17712
	
	## model on outcome (incl mediator)
	# first run with time-dependent vars
	a <- aalen(Surv(t, out) ~ exposed + carstairs + factor(cal_period) + sleep,
						 id = df_test$patid,
						 #clusters = df_test$setid,
						 data = df_test,
						 robust = T)
	summary(a) ## none of the p-values show evidence against H_0 of constant effect
	##
	b <- aalen(Surv(t, out) ~ const(exposed) + const(carstairs) + const(sleep),
						 id = df_test$patid,
						 #clusters = df_test$setid,
						 data = df_test,
						 robust = T)
	summary(b) 
	##
	c <- aalen(Surv(t, out) ~ const(exposed) + const(carstairs),
						 id = df_test$patid,
						 #clusters = df_test$setid,
						 data = df_test,
						 robust = T)
	summary(c)
	##
	full_model_cov <- b$var.gamma
	direct_effect <- tibble(names = rownames(b$gamma),
													est = as.numeric(b$gamma) * 1e4) 
	# Direct effect: exposure increases depression by 1.6 per 10,000 p-years
	
	## calculate total/direct/indirect effect
		## plus CI
	CI_comp <- function(mean_lambda1 ,
											mean_lambda3,
											covar11,
											covar12,
											covar22,
											mean_alpha ,
											var_alpha,
											G = 10 ^ 4)
	{
		require(mvtnorm)
		Omega <- matrix(c(covar11, covar12, covar12, covar22), nrow = 2)
		IE <- rep(0, G)
		TE <- rep(0, G)
		Q <- rep(0, G)
		for (i in 1:G)
		{
			lambda <- rmvnorm(1,
												mean = c(mean_lambda1, mean_lambda3),
												sigma = Omega)
			alpha <- rnorm(1, mean = mean_alpha, sd = sqrt(var_alpha))
			IE[i] <- lambda[2] * alpha
			TE[i] <- IE[i] + lambda[1]
			Q[i] <- IE[i] / TE[i]
		}
		print("IE:")
		print(mean(IE))
		print(quantile(IE, c(0.025, 0.975)))
		print("TE:")
		print(mean(TE))
		print(quantile(TE, c(0.025, 0.975)))
		print("Q:")
		print(mean(Q))
		print(quantile(Q, c(0.025, 0.975)))
		
		IE_stats <- c(mean(IE), quantile(IE, c(0.025, 0.975)))
		TE_stats <- c(mean(TE), quantile(TE, c(0.025, 0.975)))
		Q_stats <- c(mean(Q), quantile(Q, c(0.025, 0.975)))
		tibble(
			IE = IE_stats[1],
			IE_lci = IE_stats[2],
			IE_uci = IE_stats[3],
			TE = TE_stats[1],
			TE_lci = TE_stats[2],
			TE_uci = TE_stats[3],
			Q = Q_stats[1],
			Q_lci = Q_stats[2],
			Q_uci = Q_stats[3],
		)
	}
	## test example from paper
	CI_comp(mean_lambda1=0.000561 ,mean_lambda3=0.000234,
				covar11=0.000197^2, covar12=-1.12*10^(-9), covar22=0.000054^2,
				mean_alpha=0.67, var_alpha=0.066^2)

	sleep_cov_index <- which(rownames(b$gamma)=="const(sleep)")
	mediation_analysis <- CI_comp(mean_lambda1=b$gamma[1] ,mean_lambda3=b$gamma[sleep_cov_index],
					covar11=(b$var.gamma)[1,1], covar12=(b$var.gamma)[1,sleep_cov_index], covar22=(b$var.gamma)[sleep_cov_index,sleep_cov_index],
					mean_alpha=test$coefficients[2], var_alpha=vcov(test)[2,2])
	mediation_analysis <- mediation_analysis %>%
		mutate_at(vars(1:6), ~.*10^4) 
	
	
	## OR using the aalen model for the exposure~mediator + conf model:
	mediation_analysis_v2 <- CI_comp(mean_lambda1=b$gamma[1] ,mean_lambda3=b$gamma[sleep_cov_index],
																covar11=(b$var.gamma)[1,1], covar12=(b$var.gamma)[1,sleep_cov_index], covar22=(b$var.gamma)[sleep_cov_index,sleep_cov_index],
																mean_alpha=mediator_model$gamma[1], var_alpha=mediator_model$var.gamma[1,1])
	mediation_analysis_v2 <- mediation_analysis_v2 %>%
		mutate_at(vars(1:6), ~.*10^4) 
	
	
	cox
	summary(test)
	summary(b)
	mediation_analysis
	mediation_analysis_v2
	
	# Psoriasis ~ depression example: 
	# What is going on? Why is sleep less likely in exposed under additive hazard model?
	# Is carstairs/calendar period doing "too much" explanation of sleep? 
	# Especially since logistic regression says exposed incraeses sleep! Could be fucked because of data set up though
	# was something weird happening with sleep 2008-14? Seems to attenuate a lot

#}