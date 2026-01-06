data alzheimer25;
    set "~/Project 3 LDA/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

proc mi data=alzheimer25 seed=1234 out=alzheimer25_mi simple nimpute=10 round=1;
var age bmi sex adl taupet0 taupet1 taupet2 taupet3 taupet4 taupet5 taupet6
abpet0 abpet1 abpet2 abpet3 abpet4 abpet5 abpet6
cdrsb0 cdrsb1 cdrsb2 cdrsb3 cdrsb4 cdrsb5 cdrsb6;
/*by patid;*/
run;

data alzheimer_long;
    set alzheimer25_mi; /* Prendi il dataset originale (WIDE) */
    
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

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
    AGE_STD = (AGE - m_age) / s_age;
    BMI_STD = (BMI - m_bmi) / s_bmi;
run;

ods trace on;

/* 1. Esegui la GENMOD (come fatto prima) */
proc genmod data = alzheimer_long_centered descending;
    class PATID SEX(ref = '1') TIMECLSS / param = ref;
    by _imputation_;
    model CDRSB_CAT = TIME TIME * TIME SEX BMI_STD ADL 
    TIME * BMI_STD TIME * ADL 
    TIME * TIME * ADL / dist=binomial link=logit;
    repeated subject=PATID / withinsubject=TIMECLSS type=UN modelse COVB;
    ods output GEEEmpPEst = gmparms 
               ParmInfo   = gmpinfo 
               GEERCov    = gmcovb;
run;

/* 2. Sincronizzazione dei parametri per MIANALYZE */
data gmpinfo_clean;
    set gmpinfo;
    /* Rimuove l'intercetta se è quella a causare il disallineamento */
    if parameter = 'Prm5' then delete; 
    /* Rinumeriamo per far combaciare Prm1 della covarianza con il primo effetto */
run;

/* 3. Pooling senza Intercept */
proc mianalyze parms=gmparms covb=gmcovb parminfo=gmpinfo wcov bcov tcov;
    modeleffects Intercept TIME TIME * TIME SEX BMI_STD ADL 
    TIME * BMI_STD TIME * ADL 
    TIME * TIME * ADL ;
    ods output ParameterEstimates=mianalyze_results;
run;


/* Trasponiamo ignorando i nomi originali per evitare conflitti con gli asterischi */
proc transpose data=mianalyze_results out=final_betas(drop=_name_ _label_) prefix=Beta;
    var Estimate;
run;

data plot_mianalyze;
    if _n_ = 1 then set final_betas;
    
    /* Variabili di riferimento */
    SEX_val = 0; 
    BMI_val = 0; 
    ADL_val = 10;

    do TIME = 0 to 6 by 0.1;
        /* Usiamo i Beta in ordine di apparizione nel modello */
        eta = Beta1                   /* Intercept */
            + (Beta2 * TIME)          /* TIME */
            + (Beta3 * TIME**2)       /* TIME*TIME */
            + (Beta4 * SEX_val)       /* SEX */
            + (Beta5 * BMI_val)       /* BMI_STD */
            + (Beta6 * ADL_val)       /* ADL */
            + (Beta7 * TIME * BMI_val) /* TIME*BMI_STD */
            + (Beta8 * TIME * ADL_val) /* TIME*ADL */
            + (Beta9 * TIME**2 * ADL_val); /* TIME*TIME*ADL */

        Pred_Prob = exp(eta) / (1 + exp(eta));
        output;
    end;
run;

proc sgplot data=plot_mianalyze;
    title "Rischio di Demenza Grave - Risultati Combinati (MI)";
    series x=TIME y=Pred_Prob / lineattrs=(color=red thickness=2);
    yaxis label="Probabilità Predetta" grid min=0 max=1;
    xaxis label="Anni" grid;
run;



/* A. CALCOLO DELLA PROPORZIONE OSSERVATA (Dati Reali) */
proc means data=alzheimer_long nway noprint;
    class TIME;
    var CDRSB_CAT;
    output out=observed_data mean=Obs_Prop;
run;

/* B. UNIONE DEI DATASET */
/* Assicurati che nel DATA plot_mianalyze la variabile si chiami Predicted_Prob */
data plot_final_combined;
    set plot_mianalyze(in=a) observed_data(in=b);
    
    /* Etichetta per la legenda (opzionale ma utile) */
    if a then Tipo = "Predetto dal Modello";
    else if b then Tipo = "Osservato nei Dati";
run;

/* C. IL PLOT DEFINITIVO */
proc sgplot data=plot_final_combined;
    title "Validazione Modello MI-GEE: Osservato vs Predetto";
    
    /* 1. La linea dei dati reali (OSSERVATI) */
    /* Usiamo i marker (punti) per i dati reali */
    scatter x=TIME y=Obs_Prop / markerattrs=(symbol=circlefilled color=black size=10) 
            legendlabel="Proporzione Osservata";
    
    /* 2. La curva del modello (PREDETTI) */
    /* Usiamo una linea continua per la predizione teorica */
    series x=TIME y=Pred_Prob / lineattrs=(color=red thickness=3) 
           legendlabel="Probabilità Predetta (MI-GEE)";

    yaxis label="Probabilità di Demenza Grave (CDRSB >= 10)" grid min=0.1 max=0.8;
    xaxis label="Anni di Follow-up" grid values=(0 to 6 by 1);
run;


/* 1. Ordina entrambi i dataset per TIME per sicurezza */
proc sort data=plot_mianalyze; by TIME; run;
proc sort data=observed_data; by TIME; run;

/* 2. Unisci i dati (MERGE) */
data plot_final_combined;
    merge plot_mianalyze observed_data;
    by TIME;
run;

/* 3. Plot corretto */
proc sgplot data=plot_final_combined;
    title "Validazione Modello MI-GEE: Osservato vs Predetto";
    
    /* Marker per i dati reali - appariranno solo dove TIME è intero (0,1,2...) */
    scatter x=TIME y=Obs_Prop / markerattrs=(symbol=circlefilled color=black size=10) 
            legendlabel="Proporzione Osservata";
    
    /* Linea continua per la predizione - fluida grazie ai passi da 0.1 */
    series x=TIME y=Pred_Prob / lineattrs=(color=red thickness=3) 
            legendlabel="Probabilità Predetta (MI-GEE)";

    yaxis label="Probabilità di Demenza Grave (CDRSB >= 10)" grid min=0 max=1;
    xaxis label="Anni di Follow-up" grid values=(0 to 6 by 1);
run;
/* COMMENTO: Non so se sia giusta la tecnica
Però comunque almeno questo va come primo inizio */

/* TO DO: Apply this to GLMM procedure
(should be more or less the same just change the code) */

/* Another idea could be to study the whole model not reduced and then the reduced model */

/* Plot */

/* 1. Prepariamo il dataset originale in formato LONG */



data alzheimer_long_ORIG;
    set alzheimer25; /* Il dataset caricato all'inizio, prima della PROC MI */
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    do TIME = 0 to 6;
        CDRSB = cdrsb_arr[TIME];
        if not missing(CDRSB) then do;
            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;
        end;
        else CDRSB_CAT = .; 
        output;
    end;
    keep PATID TIME CDRSB_CAT;
run;

/* 2. Calcoliamo le proporzioni sui dati ORIGINALI */
proc means data=alzheimer_long_ORIG nway noprint;
    class TIME;
    var CDRSB_CAT;
    output out=stats_orig mean=Prop_Originale;
run;

/* 3. Calcoliamo le proporzioni sui dati IMPUTATI (già presente nel tuo codice, lo rinominiamo) */
proc means data=alzheimer_long nway noprint;
    class TIME;
    var CDRSB_CAT;
    output out=stats_imput mean=Prop_Imputata;
run;

/* 4. Unione dei dati per il plot */
data plot_confronto;
    merge stats_orig(keep=TIME Prop_Originale) 
          stats_imput(keep=TIME Prop_Imputata)
          plot_mianalyze(keep=TIME Pred_Prob); /* La linea del modello */
    by TIME;
run;

/* 5. Plot di confronto */
proc sgplot data=plot_confronto;
    title "Confronto Proporzioni: Dati Originali vs Imputati vs Modello";
    
    /* Proporzione sui dati reali (punti blu) */
    scatter x=TIME y=Prop_Originale / markerattrs=(symbol=circlefilled color=blue size=10) 
            legendlabel="Osservato (Dati Originali)";
            
    /* Proporzione sui dati imputati (punti neri vuoti) */
    scatter x=TIME y=Prop_Imputata / markerattrs=(symbol=circle color=black size=12) 
            legendlabel="Medie post-Imputazione (MI)";
            
    /* Linea del modello (rossa) */
    series x=TIME y=Pred_Prob / lineattrs=(color=red thickness=3) 
            legendlabel="Predizione Modello MI-GEE";

    yaxis label="Probabilità Demenza Grave" grid min=0 max=1;
    xaxis label="Anni" grid values=(0 to 6 by 1);
run;
