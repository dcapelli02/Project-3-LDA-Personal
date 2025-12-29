data alzheimer25;
    set "~/Project 3 LDA/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

data alzheimer_long;
    set alzheimer25; /* Prendi il dataset originale (WIDE) */
    
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

/* CONTROLLO DI SICUREZZA */
proc freq data=alzheimer_long_centered;
    tables R;
run;

/*%dropout(data=alzheimer_long, id=patid, time=time, response=cdrsb_cat, out=alz_dropout);*/


/* 1. Crea l'indicatore di risposta R */
data weights_prep;
    set alzheimer_long_centered;
    by patid;
    R = (CDRSB_CAT ne .);
    /* Recupera il valore della risposta alla visita precedente */
    prev_resp = lag(CDRSB_CAT);
    if first.patid then prev_resp = 0; /* o un valore di baseline */
run;

/* 2. Stima la probabilità di essere osservati (Missingness Model) */
proc logistic data=weights_prep descending;
    class SEX TIMECLSS;
    model R = TIME SEX BMI prev_resp / link=logit;
    output out=weights_data p=prob_obs;
run;

/* 3. Calcola il peso (Inverse Probability Weight) */
data alzheimer_weighted;
    set weights_data;
    if R = 1 then ipw = 1 / prob_obs;
    else delete; /* Le WGEE usano solo i casi osservati ma pesati */
run;

proc means data=alzheimer_weighted n nmiss min;
    var CDRSB_CAT TIME SEX BMI_STD ADL ipw;
run;

/* 4. Esegui il codice */
proc genmod data = alzheimer_weighted descending;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT = TIME TIME*TIME SEX BMI_STD ADL 
                      TIME*BMI_STD TIME*ADL 
                      TIME*TIME*ADL /
                      dist=binomial link=logit;

    /* L'istruzione chiave per trasformarla in WGEE */
    weight ipw;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
    output out = gee_preds p = Predicted_Prob;
run;
