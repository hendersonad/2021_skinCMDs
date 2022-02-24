library(data.table)
depression_static <- readRDS(paste0(datapath, "out/", ABBRVexp, "-depression_static.rds"))

# Who gets depression?  ---------------------------------------------------
dt_dep = data.table(depression_static)
dt_dep[, age_x := as.numeric(age_cat)]
ages = levels(dt_dep$age_cat)

dt_dep[, Age2 := as.character(age_cat)]
ages2 = unique(dt_dep[order(age)]$Age2)
dt_dep[, Age2 := factor(Age2, levels = ages2)]
dt_dep[, age2_x := as.numeric(Age2)]

dt_dep[, time := as.numeric(tstop-tstart)]
dt_dep[, dep := out]

d_standard = dt_dep[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x)]
d_standard[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
d_standard[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]

# Deaths by "what"
dc = dt_dep[!is.na(exposed)]


dc[, what := gender]#factor(ifelse(gender==1, "Male", ifelse(gender==2, "Female", "Other/not-specified")), levels = c("Male", "Female", "Other/not-specified"))]
dc[, what2 := bmi_cat]
dc[, what3 := as.factor(eth_edited)]
dc[, what4 := as.factor(alc)]

dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, ABBRVexp, "Control"), levels = c(ABBRVexp, "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

## WHAT
pdf(file = here::here("out", paste0("desc_",ABBRVexp,"_dep_gender.pdf")), 6,6)
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
	scale_shape_manual(values = c(ABBRVexp = 5, "Control" = 20)) +
	labs(x =  "Age", 
			 y = "Depression per 100\nyrs of followup" , 
			 colour = "gender", shape = NULL) +
	theme_bw() +
	theme(legend.position = "top", legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
				legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
				legend.key.size = unit(0.225, "cm")) + 
	guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
				 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
				 linetype = "none")
dev.off()

## WHAT2
dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what2, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, ABBRVexp, "Control"), levels = c(ABBRVexp, "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

pdf(file = here::here("out", paste0("desc_",ABBRVexp,"_dep_BMI.pdf")), 6,6)
ggplot() +
	geom_rect(data = d_standard, aes(xmin = age2_x - 0.5, xmax = age2_x + 0.5, 
																	 ymin = rate_lo_frequentist, ymax = rate_hi_frequentist), fill = "#cccccc") +
	geom_segment(data = d_standard, aes(x = age2_x - 0.5, xend = age2_x + 0.5, 
																			y = rate, yend = rate), size = 0.25, linetype = "22") +
	geom_pointrange(data = dp, aes(x = age2_x, ymin = rate_lo_frequentist, y = rate, ymax = rate_hi_frequentist, 
																 colour = what2, shape = z_exp_label, linetype = z_exp_label), size = 1.4, fatten = 2.4, position = position_dodge(width = 1), stroke = 0.5) +
	scale_x_continuous(breaks = 1:length(ages2), labels = ages2, expand = expansion(0.01)) +
	scale_linetype_manual(values = ltys) +
	scale_colour_manual(values = clrs) +
	scale_shape_manual(values = c(ABBRVexp = 5, "Control" = 20)) +
	labs(x =  "Age", 
			 y = "Depression per 100\nyrs of followup" , 
			 colour = "BMI", shape = NULL) +
	theme_bw() +
	theme(legend.position = "top", legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
				legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
				legend.key.size = unit(0.225, "cm")) + 
	guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
				 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
				 linetype = "none")
dev.off()

## WHAT3
dp = dc[, .(dep = sum(dep), years = sum(time), rate = 100 * sum(dep) / sum(time)), keyby = .(age2_x, what3, exposed)]
dp[, rate_lo_frequentist := qchisq(0.025, 2 * dep) / (2 * years / 100)]
dp[, rate_hi_frequentist := qchisq(0.975, 2 * dep + 2) / (2 * years / 100)]
dp[, z_exp_label := factor(ifelse(exposed == 1, ABBRVexp, "Control"), levels = c(ABBRVexp, "Control"))]

clrs = c("#6388b4", "#fd57ca", "#47a74e", "#eb1e2c", "#9c5142", "#8175aa", "#ccb22b")
ltys = 1:5

pdf(file = here::here("out", paste0("desc_",ABBRVexp,"_dep_eth.pdf")), 6,6)
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
	scale_shape_manual(values = c(ABBRVexp = 5, "Control" = 20)) +
	labs(x =  "Age", 
			 y = "Depression per 100\nyrs of followup" , 
			 colour = "Eth5", shape = NULL) +
	theme_bw() +
	theme(legend.position = "top", legend.title = element_text(size = 9, margin = margin(0, 0, -0.05, 0, "cm")), 
				legend.text = element_text(size = 9, margin = margin()), legend.justification = c(0, 1),
				legend.key.size = unit(0.225, "cm")) + 
	guides(colour = guide_legend(order = 1, override.aes = list(size = 0.4)), 
				 shape = guide_legend(order = 2, reverse = TRUE, override.aes = list(size = 0.4)),
				 linetype = "none")
dev.off()

