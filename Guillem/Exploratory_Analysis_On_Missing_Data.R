library(haven)      # Dades SAS
library(dplyr)      # Manipulació Dades
library(tidyr)      # Transformacions (pivot_longer)
library(ggplot2)    # Gràfics
library(naniar)     # Missing Data 
library(VIM)        # Visualització Patrons
library(nlme)       # LMM 
library(lme4)       # GLMM 
library(geepack)    # GEE i WGEE 
library(mice)       # Sensitivity Analysis

#================================================================#

# DATASET

alz <- read_sas("C:/Users/win11/Documents/Guillem/Erasmus/Assignatures/Longitudinal Data Analysis/Project-3-LDA-Personal/alzheimer25.sas7bdat")

# Factor variables

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

# Bins

alz$cdrsb_bin0 <- ifelse(alz$cdrsb0 > 10, 1, 0)
alz$cdrsb_bin1 <- ifelse(alz$cdrsb1 > 10, 1, 0)
alz$cdrsb_bin2 <- ifelse(alz$cdrsb2 > 10, 1, 0)
alz$cdrsb_bin3 <- ifelse(alz$cdrsb3 > 10, 1, 0)
alz$cdrsb_bin4 <- ifelse(alz$cdrsb4 > 10, 1, 0)
alz$cdrsb_bin5 <- ifelse(alz$cdrsb5 > 10, 1, 0)
alz$cdrsb_bin6 <- ifelse(alz$cdrsb6 > 10, 1, 0)

# Baseline

alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0
alz$bprs_base <- alz$bprs0

#================================================================#

# LONGITUDINAL DATASET

alz_long <- alz %>%
  pivot_longer(
    cols = matches(".*\\d+$"),     # Selecciona qualsevol columna que acabi en número
    names_to = c(".value", "time"),
    names_pattern = "(.*)(\\d+)$"  
  ) %>%
  mutate(
    time = as.numeric(time),       
    id = as.factor(patid)         
  ) %>%
  arrange(id, time)

head(alz_long)
str(alz_long)

#================================================================#

# EXPLORATORY ANALYSIS ON MISSING DATA

# General Behavior

alz %>%
  select(matches("\\d+$")) %>% 
  vis_miss() +
  labs(title = "Missingness Pattern") +
  theme(axis.text.x = element_text(angle = 90))

alz %>%
  select(matches("cdrsb\\d+$")) %>% 
  vis_miss() +
  labs(title = "Missingness Pattern on cdrsb") +
  theme(axis.text.x = element_text(angle = 90))

alz_long %>%
  select(cdrsb, bprs, abpet, taupet) %>% 
  gg_miss_var(show_pct = TRUE) + 
  labs(title = "% Missing Data",
       y = "Missings' Percentage")

cat("The missingness pattern is equal in every longitudinal variable.")

# Dropout

dropout_table <- alz_long %>%
  group_by(time) %>%
  summarise(
    N_Total = n(),
    Observed_CDRSB = sum(!is.na(cdrsb)),
    Missing_CDRSB = sum(is.na(cdrsb)),
    Pct_Missing = round(mean(is.na(cdrsb)) * 100, 2)
  )

print(dropout_table)

ggplot(dropout_table, aes(x = time, y = Observed_CDRSB)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Number of observed patients evolution",
       x = "Year", y = "Observed Patients") +
  theme_bw()

# Monotoneness

check_intermittent <- alz_long %>%
  select(id, time, cdrsb) %>%
  arrange(id, time) %>%
  group_by(id) %>%
  summarise(
    is_missing = list(is.na(cdrsb)),
    pattern = paste(as.integer(!is.na(cdrsb)), collapse = "")
  )
table(check_intermittent$pattern)

cat("We conclude that missingness follows a monotone pattern.")

# Missingness per group

# Trial
alz %>%
  select(trial, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = trial) +
  labs(title = "% Missings depending on Trial",
       y = "% Observed Patients",
       x = "Trial")

# Sex
alz %>%
  select(sex, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = sex) +
  labs(title = "% Missings depending on Sex",
       y = "% Observed Patients",
       x = "Sex")

# Age
alz %>%
  select(age, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = age) +
  labs(title = "% Missings depending on Age",
       y = "% Observed Patients",
       x = "Age")

# Edu
alz %>%
  select(edu, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = edu) +
  labs(title = "% Missings depending on Edu",
       y = "% Observed Patients",
       x = "Edu")

# BMI
alz %>%
  select(bmi, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = bmi) +
  labs(title = "% Missings depending on BMI",
       y = "% Observed Patients",
       x = "BMI")

# Inkomen
alz %>%
  select(inkomen, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = inkomen) +
  labs(title = "% Missings depending on Inkomen",
       y = "% Observed Patients",
       x = "Inkomen")

# Job
alz %>%
  select(job, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = job) +
  labs(title = "% Missings depending on Job",
       y = "% Observed Patients",
       x = "Job")

# ADL
alz %>%
  select(adl, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = adl) +
  labs(title = "% Missings depending on ADL",
       y = "% Observed Patients",
       x = "ADL")

# WZC
alz %>%
  select(wzc, matches("cdrsb\\d+$")) %>%  
  gg_miss_fct(fct = wzc) +
  labs(title = "% Missings depending on WZC",
       y = "% Observed Patients",
       x = "WZC")


