

# testing interactions  ---------------------------------------------------
exposure <- "psoriasis"
ABBRVexp <- substr(exposure, 1, 3)
outcome <- "depression"

df_model <- readRDS(paste0(datapath,"out/models_data/df_model",ABBRVexp,"_",outcome,".rds"))

mod5 <- readRDS(paste0(datapath, "out/models_data/", ABBRVexp, "_", outcome, "_mod5_interaction_age_modeldata.rds"))

mod5$call

xx_ids <- unique(df_model$patid)

set.seed(1435)
xx_sample <- sample(xx_ids, size = 10000)
df_sample <- df_model %>% 
  filter(setid %in% xx_sample)

levels(df_sample$agegroup)
cox1 <- coxph(formula = Surv(t, out) ~ exposed + agegroup + exposed * 
                agegroup + carstairs + cal_period + comorbid + cci + strata(setid), 
              data = df_sample)
broom::tidy(cox1, conf.int = TRUE, conf.level = 0.99, exp = T) %>% print(n=Inf)

old_lvl <- levels(df_sample$agegroup)
new_lvl <- old_lvl[c(2,3,1,4,5,6)]

df_sample$agegroup <- as.character(df_sample$agegroup)
df_sample$agegroup <- factor(df_sample$agegroup, levels = new_lvl)
cox2 <- coxph(formula = Surv(t, out) ~ exposed + agegroup + exposed * 
                agegroup + carstairs + cal_period + comorbid + cci + strata(setid), 
              data = df_sample)


print("Estimating HR in youngest group from model 1 (agegroup50-59 as ref group) - 2 options")
lincom(cox1, 
       c("exposedPsoriasis+exposedPsoriasis:agegroup18-39",
         "exposedPsoriasis+agegroup18-39+exposedPsoriasis:agegroup18-39"),
       level = 0.99, eform = T)

print("Switch ref group to the youngest and print the maine exposure HR")
broom::tidy(cox2, conf.int = TRUE, conf.level = 0.99, exp = T) %>% filter(term == "exposedPsoriasis") 

print("Matches the exposedPsoriasis+exposedPsoriasis:agegroup18-39 paramaterisation (i.e. without the agegroup coefficient)")


# testing predicted vals stuff --------------------------------------------
predicted_vals <- predict(mod5, type = "expected", se.fit = F)
predicted_vals <- as.vector(predicted_vals)

dim(df_model)
df_cox <- df_model %>% 
  dplyr::select(setid, patid, t, out, exposed, agegroup, carstairs, cal_period, comorbid, cci) %>% 
  drop_na()
dim(df_cox)
df_cox$pred <- predicted_vals
df_pred <- df_cox %>% 
  left_join(dplyr::select(df_sample, setid, patid, age)) %>% 
  filter(pred > 0.001)

ggplot(df_pred, aes(x = pred, group = exposed, colour = exposed, fill=exposed)) +
  geom_density(alpha = 0.2) + 
  theme_ali()

ggplot(df_pred, aes(x = t, y = pred, group = exposed, colour = exposed, fill=exposed)) +
  geom_point(alpha = 0.2) + 
  theme_ali()
