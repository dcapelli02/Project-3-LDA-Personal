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

/* LME - Linear Mixed-Effects Model */

proc mixed data=alzheimer_long_centered method=REML covtest plots(maxpoints=none);
    class PATID TRIAL SEX EDU JOB WZC / ref=first; 
    
    model BPRS = TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
    			TIME * AGE TIME * EDU TIME * JOB TIME * ADL TIME * CDRSB0/ solution cl; 
    

    random intercept TIME / subject=PATID type=UN g;
    
    title "Model 1: Linear Mixed-Effects Model";
run;

/* ============================================================== */

/* Sensitivity Analysis to access MNAR */


/* MNAR */

/* Multiple Imputation (MI) - M = 10 */

proc mi data=alzheimer25 out=mi_mar_temp nimpute=10 seed=11;
    class SEX TRIAL JOB WZC EDU; 
    
    var BPRS0 BPRS1 BPRS2 BPRS3 BPRS4 BPRS5 BPRS6 
        AGE BMI EDU JOB ADL WZC CDRSB0 SEX TRIAL;
        
    fcs; 
    title "Sensitivity Analysis: Step 1 - Standard Imputation";
run;



proc sort data=alzheimer25; by PATID; run;
proc sort data=mi_mar_temp; by PATID; run;

data mi_mnar_out;
    merge mi_mar_temp(in=a) 
          alzheimer25(keep=PATID BPRS1 BPRS2 BPRS3 BPRS4 BPRS5 BPRS6 
                      rename=(BPRS1=o_B1 BPRS2=o_B2 BPRS3=o_B3 
                              BPRS4=o_B4 BPRS5=o_B5 BPRS6=o_B6));
    by PATID;
    if a; 


    array imp[*] BPRS1-BPRS6;      /* Imputed Values */
    array orig[*] o_B1-o_B6;       /* Original Values (with missings) */

    do i = 1 to dim(imp);
    	/* We sum +5 to BPRS when they dropout (Shift=5) */
		/* This simulates that patients who leave are worse (higher BPRS) */
        if missing(orig[i]) then imp[i] = imp[i] + 5;
    end;
    
    drop i o_B1-o_B6; 
run;

/* Dataset in longitudinal format */
data alzheimer_mnar_long;
    if _n_ = 1 then set stats;    
    set mi_mnar_out; /* Imputed Data */
    
    array bprs_arr[0:6] bprs0-bprs6; 
    array cdrsb_arr[0:6] cdrsb0-cdrsb6; 
    
    do TIME = 0 to 6;
        BPRS = bprs_arr[TIME];
        
        AGE_STD = (AGE - mean_age) / sd_age;
        BMI_STD = (BMI - mean_bmi) / sd_bmi;
        
        TIMECLSS = put(TIME, 1.);
        
        if not missing(BPRS) then output;
    end;
    
    keep PATID _Imputation_ TIME TIMECLSS BPRS CDRSB0 
         SEX AGE AGE_STD EDU BMI BMI_STD JOB ADL WZC TRIAL;
run;

/* LME Model */

proc sort data=alzheimer_mnar_long;
    by _Imputation_ PATID; /* Ordered */
run;

proc mixed data=alzheimer_mnar_long method=REML maxiter=200 scoring=5 plots(maxpoints=none) ;
    by _Imputation_; /* This analyses the 10 datasets separately */
    class PATID TRIAL SEX EDU JOB WZC / ref=first; 
    
    model BPRS = TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
                 TIME * AGE TIME * EDU TIME * JOB TIME * ADL TIME * CDRSB0 / solution; 
    
    random intercept TIME / subject=PATID type=UN;
    
    ods output SolutionF=mix_parms;
    title "Fitting LME Model on MNAR Imputed Data";
run;

/* IMPORTANT: Table to Compare */

proc mianalyze parms=mix_parms;
    class TRIAL SEX EDU JOB WZC;
    modeleffects Intercept TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
                 TIME*AGE TIME*EDU TIME*JOB TIME*ADL TIME*CDRSB0;
    title "Final Sensitivity Results (MNAR Shift=5)";
run;

/* We compare the p-values of the normal LME and the LME with Sensitivity Analysis.
We look if they are similar or not: if they are similar, it means that the results are robust;
if they are not similar, it means that the results are sensitive.
In our case: TIME * CDRSB0, TIME * AGE, TIME * EDU, AGE, BMI and WZC are robust
and TIME * JOB, TIME * ADL and CDRSB0 are sensitive. 
Although the initial analysis suggested that JOB and ADL influenced disease progression, 
the sensitivity analysis indicates that these associations may be biased by patient dropout. */
