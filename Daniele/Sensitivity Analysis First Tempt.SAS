data alzheimer25;
    set "~/Project 3 LDA/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;


/* ==========================================================
   STEP 1: MULTIPLE IMPUTATION (MNAR SENSITIVITY - SHIFT METHOD)
   ========================================================== */
/* Ipotizziamo uno scenario MNAR in cui chi abbandona ha un 
   peggioramento di +2 punti nel CDRSB rispetto al previsto */

proc mi data=alzheimer25 seed=1234 out=alzheimer_mi_mnar nimpute=10 round=1;
    var age bmi sex bprs0 bprs1 bprs2 bprs3 bprs4 bprs5 bprs6;
    fcs reg; 
    /* Applichiamo lo shift per la Sensitivity Analysis */
    mnar adjust (bprs1 / shift=5)
         adjust (bprs2 / shift=5)
         adjust (bprs3 / shift=5)
         adjust (bprs4 / shift=5)
         adjust (bprs5 / shift=5)
         adjust (bprs6 / shift=5);
run;

/* ==========================================================
   STEP 2: TRASFORMAZIONE DA WIDE A LONG
   ========================================================== */
data alzheimer_long;
    set alzheimer_mi_mnar; /* Prendi il dataset originale (WIDE) */
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];
        
        /* 1. CREA L'INDICATORE DI RISPOSTA */
        /* Se CDRSB è vuoto (.), R sarà 0. Altrimenti 1. */
        if missing(CDRSB) then R = 0; 
        else R = 1;

        /* 2. CREA LA CATEGORIA SOLO SE IL DATO ESISTE */
        if R = 1 then do;
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;
        end;
        else CDRSB_CAT = .; /* Importante: lascia vuoto se manca il dato */

        TIMECLSS = put(TIME, 1.);
        
        /* 3. L'ISTRUZIONE CRUCIALE */
        /* Questo 'output' deve essere FUORI da ogni 'if'. 
           Deve scrivere una riga per ogni iterazione del loop 'do'. */
        output; 
    end;
    
    drop cdrsb0 cdrsb1 cdrsb2 cdrsb3 cdrsb4 cdrsb5 cdrsb6
    abpet0 abpet1 abpet2 abpet3 abpet4 abpet5 abpet6
    taupet0 taupet1 taupet2 taupet3 taupet4 taupet5 taupet6;
run;

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET AGE BMI; 
    output out=bio_stats mean=m_ab m_tau m_age m_bmi std=s_ab s_tau s_age s_bmi;
run;


/* ==========================================================
   STEP 3: STANDARDIZZAZIONE E CENTRATURA (Sui dati imputati)
   ========================================================== */
proc sort data=alzheimer_long; by _Imputation_; run;

proc means data=alzheimer_long noprint;
    by _Imputation_;
    var age bmi;
    output out=stats_mi mean=m_age m_bmi std=s_age s_bmi;
run;

data alzheimer_long_std;
    merge alzheimer_long stats_mi;
    by _Imputation_;
    age_std = (age - m_age) / s_age;
    bmi_std = (bmi - m_bmi) / s_bmi;
run;

/* ==========================================================
   STEP 4: ANALISI DEI MODELLI MISTI (BY IMPUTATION)
   ========================================================== */
/* Inseriamo anche l'analisi di INFLUENCE come extra */

ods exclude all;

proc mixed data=alzheimer_long_std method=REML;
    by _Imputation_;
    class sex patid;
    model bprs = age_std bmi_std sex time age_std*time 
                  / solution influence(iter=5);
    random intercept time / subject=patid type=un;
    ods output SolutionF=estimates_mnar;
run;

ods select all;

/* ==========================================================
   STEP 5: POOLING CON MIANALYZE (RISULTATO FINALE)
   ========================================================== */
proc mianalyze parms=estimates_mnar;
    modeleffects Intercept age_std bmi_std sex time age_std*time;
run;
