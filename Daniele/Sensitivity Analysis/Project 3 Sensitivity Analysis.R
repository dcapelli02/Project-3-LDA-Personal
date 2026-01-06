## Libraries
library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(nlme)
library(lmtest)
library(DescTools)
library(Matrix)
library(MASS)
library(metafor)
library(geepack)
library(car)
library(lme4)
library(ORTH.Ord)
library(MuMIn)
library(patchwork)
library(naniar)
library(mice)
library(tidyverse)
library(lmerTest) # Per i p-value nei modelli misti
library(broom.mixed)
library(purr)

## Import and fix the data
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

head(alz)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))
alz$cdrsb_base <- alz$cdrsb0


run_shift_analysis <- function(data_wide, shift_value) {
  # 1. Configurazione MICE
  ini <- mice(data_wide, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  meth[paste0("bprs", 1:6)] <- "norm"
  
  post <- ini$post
  for (v in paste0("bprs", 1:6)) {
    post[v] <- paste0("imp[[j]][, i] <- imp[[j]][, i] + ", shift_value)
  }
  
  # 2. Imputazione
  imp <- mice(data_wide, m = 10, method = meth, post = post, 
              seed = 1234, ridge = 0.001, printFlag = FALSE)
  
  print("--Fine imputazione--")
  
  # 3. Trasformazione Long e Standardizzazione
  alz_long <- complete(imp, action = "long") %>%
    pivot_longer(
      cols = matches("^(bprs|)[0-6]$"),
      names_to = c(".value", "TIME"),
      names_pattern = "([a-z]+)([0-6])"
    ) %>%
    mutate(year = as.numeric(TIME)) %>%
    group_by(.imp) %>%
    mutate(
      age_std = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
      bmi_std = (bmi - mean(bmi, na.rm = TRUE)) / sd(bmi, na.rm = TRUE)
    ) %>%
    ungroup()
  
  print("--Fine creazione dataset--")
  #View(alz_long)
  
  # 4. Fit del modello (usando lme4 direttamente su oggetto mira per semplicità)
  # Usiamo with() che è il comando standard di mice per modelli su dati imputati
  fit <- alz_long %>%
    group_by(.imp) %>%
    do(model = lmer(bprs ~ trial + age + edu + bmi + job + adl_num + 
                      wzc + cdrsb_base + year + age:year + edu:year + job:year + 
                      adl_num:year + cdrsb_base:year + (1 + year | .id), # Intercetta casuale
                    data = ., REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa")))
  
  # Estrai la lista pura dei modelli
  models_only <- as.list(fit$model)
  
  # Nomina i modelli nella lista (opzionale ma aiuta)
  names(models_only) <- 1:length(models_only)
  
  # Pooling
  pooled_results <- pool(as.mira(models_only))
}

shifts_da_testare <- c(0, 1, 2, 3, 4, 5)
nomi_scenari <- c("MAR (Shift 0)", "MNAR Lieve (Shift 1)",
                  "MNAR Lieve (Shift 2)",  "MNAR Tostino (Shift 3)",
                  "MNAR Fortino (Shift 4)", "MNAR Forte (Shift 5)")

# Eseguiamo le analisi
risultati_confronto <- lapply(shifts_da_testare, function(s) {
  message(paste("Eseguendo analisi con shift:", s))
  run_shift_analysis(alz, s)
})

names(risultati_confronto) <- nomi_scenari


confronto_tabella <- map_df(nomi_scenari, function(nome) {
  summary(risultati_confronto[[nome]], conf.int = TRUE) %>%
    mutate(scenario = nome)
})

tabella_stima <- confronto_tabella %>%
  dplyr::select(term, estimate, scenario) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = estimate)

# Visualizza
View(tabella_stima)


tabella_confronto_completa <- confronto_tabella %>%
  # Creiamo una stringa che contiene Stima e (std error)
  mutate(valore = sprintf("%.3f (sd=%.3f)", estimate, std.error)) %>%
  dplyr::select(term, scenario, valore) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = valore)

View(tabella_confronto_completa)



#### NCMV Pattern -- No shift ####

# 1. Identifichiamo il pattern di ogni paziente (ultimo tempo osservato)
# Assumiamo bprs0-bprs6 siano nelle colonne 18:24
alz$pattern <- apply(alz[, paste0("bprs", 0:6)], 1, function(x) max(which(!is.na(x))) - 1)
alz$pattern <- as.factor(alz$pattern) # Fondamentale che sia factor

# 2. Setup MICE
ini <- mice(alz, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# 3. Logica NCMV (Temporalità + Pattern)
vars_bprs <- paste0("bprs", 0:6)

for (i in 2:length(vars_bprs)) {
  current_var <- vars_bprs[i]
  
  # A. Solo passato (Temporalità come avevi fatto tu)
  future_vars <- vars_bprs[(i+1):length(vars_bprs)]
  future_vars <- future_vars[!is.na(future_vars)]
  if (length(future_vars) > 0) pred[current_var, future_vars] <- 0
  
  # B. Forza l'uso della variabile 'pattern' per questa visita
  # Questo dice a MICE di imputare la bprs basandosi sul momento in cui il paziente è uscito
  pred[current_var, "pattern"] <- 1 
  
  meth[current_var] <- "norm"
}

# 4. Esecuzione (Senza Shift o con Shift se vuoi testare MNAR pesante)
imp_ncmv_style <- mice(alz, 
                       m = 10, 
                       predictorMatrix = pred, 
                       method = meth, 
                       seed = 486048,
                       printFlag = FALSE)

# ==========================================================
# STEP 2: TRASFORMAZIONE DA WIDE A LONG
# ==========================================================
# Estraiamo i dati imputati (long format rispetto alle iterazioni m)
alz_long_m <- complete(imp_ncmv_style, action = "long")

# Pivot dei dati per renderli longitudinali (long format rispetto al tempo)
alz_long_final <- alz_long_m %>%
  pivot_longer(
    cols = matches("^(bprs|)[0-6]$"),
    names_to = c(".value", "TIME"),
    names_pattern = "([a-z]+)([0-6])"
  ) %>%
  mutate(
    year = as.numeric(TIME),
    # Creazione indicatori come in SAS
    TIMECLSS = as.factor(TIME)
  )

# ==========================================================
# STEP 3: STANDARDIZZAZIONE E CENTRATURA (Sui dati imputati)
# ==========================================
alz_long_std <- alz_long_final %>%
  group_by(.imp) %>%
  mutate(
    age_std = (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE),
    bmi_std = (bmi - mean(bmi, na.rm = TRUE)) / sd(bmi, na.rm = TRUE)
  ) %>%
  ungroup()

# 6. ANALISI DEI MODELLI MISTI (STEP 4)
fit_list <- alz_long_std %>%
  filter(.imp > 0) %>%
  group_by(.imp) %>%
  do(model = lmer(bprs ~ trial + age + edu + bmi + job + adl_num + 
                    wzc + cdrsb_base + year + age:year + edu:year + job:year + 
                    adl_num:year + cdrsb_base:year + (1 + year | .id), # Intercetta casuale
                  data = ., REML = TRUE,
                  control = lmerControl(optimizer = "bobyqa")))

# 7. POOLING (STEP 5)
pooled_results <- pool(as.mira(as.list(fit_list$model)))
summary(pooled_results)






# Estraiamo i risultati del primo modello e aggiungiamo l'etichetta "Modello Iniziale"
risultato_iniziale <- summary(pooled_results, conf.int = TRUE) %>%
  mutate(scenario = "NCMV")

# Uniamo questo risultato alla tabella degli shift
confronto_totale <- bind_rows(risultato_iniziale, confronto_tabella)

# Versione con Stima e Deviazione Standard (sd) come la chiamavi tu
tabella_unita_finale <- confronto_totale %>%
  mutate(valore = sprintf("%.3f (sd=%.3f)", estimate, std.error)) %>%
  dplyr::select(term, scenario, valore) %>%
  # Usiamo fct_relevel per assicurarci che il Modello Iniziale sia la prima colonna
  mutate(scenario = factor(scenario, levels = c("NCMV", nomi_scenari))) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = valore)

# Visualizza il risultato finale completo
View(tabella_unita_finale)


# Installa se non lo hai: install.packages("writexl")
library(writexl)

write_xlsx(tabella_unita_finale, "Tabella_Confronto_MNAR.xlsx")