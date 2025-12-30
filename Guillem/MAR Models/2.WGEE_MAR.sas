data alzheimer25;
    set "/home/u64400254/sasuser.v94/alzheimer25.sas7bdat"; 
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

data alzheimer_long;
    if _n_ = 1 then set stats;   
    set alzheimer25;
    
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6] bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;
    
    do TIME = 0 to 6;
        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];

        if not missing(CDRSB) 
         and not missing(BPRS)
         and not missing(ABPET)
         and not missing(TAUPET) then do;

            if CDRSB < 10 then CDRSB_CAT = 0;
            else CDRSB_CAT = 1;

            AGE_STD = (AGE - mean_age) / sd_age;
            BMI_STD = (BMI - mean_bmi) / sd_bmi;

        	TIMECLSS = put(TIME, 1.);

            output;
        end;
    end;

	keep PATID SEX AGE AGE_STD BMI BMI_STD ADL 
     	TIME TIMECLSS CDRSB CDRSB_CAT BPRS ABPET CDRSB0 TAUPET ABPET0 TAUPET0
     	TRIAL JOB WZC EDU INKOMEN;
run;

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET; 
    output out=bio_stats mean=m_ab m_tau std=s_ab s_tau;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;
    
    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
run;

/* ============================================================== */

/* GEE - Generalized Estimating Equations Model */

/* Primary Model */

proc genmod data = alzheimer_long_centered descending;
    class PATID SEX TIMECLSS;

	model CDRSB_CAT = TIME TIME * TIME SEX BMI_STD ADL 
     TIME * BMI_STD TIME * ADL  
     TIME * TIME * ADL/

        dist=binomial link=logit;

    repeated subject=PATID / withinsubject=timeclss type=UN covb corrw modelse;
    output out = gee_preds p = Predicted_Prob;
    
    title "Model 2: Generalized Estimating Equations Model";
run;

/* ============================================================== */

/* MAR Model */

proc sort data=alzheimer25; by PATID; run;

data dropout_prep;
    set alzheimer25;
    /* For each visit, we look if it's observed */
	array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    do TIME = 0 to 5;
        /* R = 1 if observed, R = 0 otherwise */
        R = (cdrsb_arr[TIME+1] ne .); 
        if cdrsb_arr[TIME] ne . then output; 
        if cdrsb_arr[TIME+1] = . then leave; 
    end;
    keep PATID TIME R SEX AGE EDU BMI WZC ADL CDRSB0;
run;

/* Dropout Model */
/* To get the weights */
proc logistic data=dropout_prep descending;
    class SEX WZC EDU;
    model R = TIME SEX AGE EDU BMI ADL WZC CDRSB0; 
    output out=probs p=prob_obs;
    title "Dropout Model";
run;

/* Shift - The prediction from this year (Y) is the one from the next year (Y+1) */ 
data probs_shifted;
    set probs;
    TIME = TIME + 1; 
    keep PATID TIME prob_obs;
run;

data wgee_data;
	merge alzheimer_long_centered (in=a) probs_shifted;
    by PATID TIME;
    if a;
    if TIME = 0 then W = 1;
    else if prob_obs > 0 then W = 1 / prob_obs;
    else W = 1; 
run;

/* WGEE */
proc genmod data=wgee_data descending;
    class PATID SEX TIMECLSS;
    model CDRSB_CAT = TIME TIME*TIME SEX BMI_STD ADL 
                      TIME*BMI_STD TIME*ADL TIME*TIME*ADL / dist=binomial link=logit;
    weight W; 
    repeated subject=PATID / withinsubject=timeclss type=UN modelse;
    title "Model 2: Weighted Generalized Estimated Equations";
run;

/* Standard GEE is only valid for MCAR, to ensure validity under MAR, we implemented WGEE 
with Inverse Probability Weighting (IPW).
The comparison between GEE and WGEE (see table) shows small differences in the estimates.
This suggest that the impact of missing data on the clinical conclusions is not as substantial as we thought.

Weights were time-shifted to ensure that the probability of dropout, 
calculated from previous visit data, was applied to the correct subsequent observation.
