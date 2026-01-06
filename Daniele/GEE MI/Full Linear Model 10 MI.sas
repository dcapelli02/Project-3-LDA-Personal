data alzheimer25;
    set "~/Project 3 LDA/alzheimer25.sas7bdat"; 
run;


/*-------------------------------------------------------------------*/
/* 1. CARICAMENTO E PREPARAZIONE DATI (WIDE TO LONG & IMPUTATION)    */
/*-------------------------------------------------------------------*/

/* Imputazione Multipla (10 set) */
proc mi data=alzheimer25 seed=1234 out=alzheimer25_mi simple nimpute=10 round=1;
    var age bmi sex adl taupet0-taupet6 abpet0-abpet6 cdrsb0-cdrsb6;
run;

/* Trasformazione in formato LONG del dataset imputato */
data alzheimer_long;
    set alzheimer25_mi;
    array cd_arr[0:6] cdrsb0-cdrsb6;
    array ab_arr[0:6] abpet0-abpet6;
    array ta_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB  = cd_arr[TIME];
        ABPET  = ab_arr[TIME];
        TAUPET = ta_arr[TIME];
        
        /* Creazione variabile target binaria */
        if not missing(CDRSB) then do;
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;
        end;
        
        TIMECLSS = put(TIME, 1.);
        output; 
    end;
    drop cdrsb0-cdrsb6 abpet0-abpet6 taupet0-taupet6;
run;

/* Standardizzazione delle covariate (su medie globali) */
proc means data=alzheimer_long noprint;
    var ABPET TAUPET AGE BMI; 
    output out=bio_stats mean=m_ab m_tau m_age m_bmi std=s_ab s_tau s_age s_bmi;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
    AGE_STD = (AGE - m_age) / s_age;
    BMI_STD = (BMI - m_bmi) / s_bmi;
run;

/*-------------------------------------------------------------------*/
/* 2. MODELLO GEE FULL LINEARE (SENZA TERMINE QUADRATICO)            */
/*-------------------------------------------------------------------*/

proc genmod data = alzheimer_long_centered descending;
    class PATID SEX(ref = '1') TIMECLSS / param = ref;
    by _imputation_;
    /* Modello Lineare: rimosso TIME*TIME */
    model CDRSB_CAT = TIME SEX AGE_STD BMI_STD ADL TAUPET_STD ABPET_STD 
                     TIME*SEX TIME*AGE_STD TIME*BMI_STD TIME*ADL 
                     TIME*TAUPET_STD TIME*ABPET_STD / 
                     dist=binomial link=logit;
    repeated subject=PATID / withinsubject=TIMECLSS type=UN modelse COVB;
    ods output GEEEmpPEst = gmparms_lin 
               ParmInfo   = gmpinfo_lin 
               GEERCov    = gmcovb_lin;
run;

/*-------------------------------------------------------------------*/
/* 3. POOLING DEI RISULTATI CON MIANALYZE                            */
/*-------------------------------------------------------------------*/

proc mianalyze parms=gmparms_lin covb=gmcovb_lin parminfo=gmpinfo_lin;
    modeleffects Intercept TIME SEX AGE_STD BMI_STD ADL TAUPET_STD ABPET_STD 
                 TIME*SEX TIME*AGE_STD TIME*BMI_STD TIME*ADL 
                 TIME*TAUPET_STD TIME*ABPET_STD;
    ods output ParameterEstimates=mianalyze_results_lin;
run;

/*-------------------------------------------------------------------*/
/* 4. COSTRUZIONE DEL PLOT DI VALIDAZIONE                            */
/*-------------------------------------------------------------------*/

/* Estrazione Beta per il calcolo della probabilità predetta */
proc transpose data=mianalyze_results_lin out=final_betas_lin(drop=_name_ _label_) prefix=Beta;
    var Estimate;
run;

data plot_mianalyze_lin;
    if _n_ = 1 then set final_betas_lin;
    /* Baseline per il plot: Maschio (0), Medie (0), ADL standard (10) */
    SEX_val = 0; AGE_val = 0; BMI_val = 0; ADL_val = 10; TAU_val = 0; ABP_val = 0;

    do TIME = 0 to 6 by 0.1;
        eta = Beta1                           /* Intercept */
            + (Beta2 * TIME)                  /* TIME */
            + (Beta3 * SEX_val)               /* SEX */
            + (Beta4 * AGE_val)               /* AGE_STD */
            + (Beta5 * BMI_val)               /* BMI_STD */
            + (Beta6 * ADL_val)               /* ADL */
            + (Beta7 * TAU_val)               /* TAUPET_STD */
            + (Beta8 * ABP_val)               /* ABPET_STD */
            + (Beta9 * TIME * SEX_val)        /* TIME*SEX */
            + (Beta10 * TIME * AGE_val)       /* TIME*AGE_STD */
            + (Beta11 * TIME * BMI_val)       /* TIME*BMI_STD */
            + (Beta12 * TIME * ADL_val)       /* TIME*ADL */
            + (Beta13 * TIME * TAU_val)       /* TIME*TAU_STD */
            + (Beta14 * TIME * ABP_val);      /* TIME*ABP_STD */

        Pred_Prob = exp(eta) / (1 + exp(eta));
        output;
    end;
run;

/* Calcolo Proporzioni Osservate sui Dati Originali (per confronto) */
data alzheimer_long_ORIG;
    set alzheimer25; 
    array cd_arr[0:6] cdrsb0-cdrsb6;
    do TIME = 0 to 6;
        CDRSB = cd_arr[TIME];
        if not missing(CDRSB) then do;
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;
        end;
        else CDRSB_CAT = .; 
        output;
    end;
run;

proc means data=alzheimer_long_ORIG nway noprint;
    class TIME;
    var CDRSB_CAT;
    output out=stats_orig mean=Prop_Originale;
run;

/* Unione e Plot Finale */
data plot_final;
    merge stats_orig(keep=TIME Prop_Originale) 
          plot_mianalyze_lin(keep=TIME Pred_Prob);
    by TIME;
run;

proc sgplot data=plot_final;
    title "Validazione: Modello Full LINEARE vs Osservato";
    scatter x=TIME y=Prop_Originale / markerattrs=(symbol=circlefilled color=blue size=10) 
            legendlabel="Osservato (Dati Originali)";
    series x=TIME y=Pred_Prob / lineattrs=(color=green thickness=3) 
            legendlabel="Predizione Modello Lineare (MI-GEE)";
    yaxis label="Probabilità Demenza Grave" grid min=0 max=1;
    xaxis label="Anni di Follow-up" grid values=(0 to 6 by 1);
run;
