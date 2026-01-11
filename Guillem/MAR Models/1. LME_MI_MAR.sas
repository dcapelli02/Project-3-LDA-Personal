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
    			TIME * AGE TIME * EDU TIME * JOB TIME * ADL TIME * CDRSB0/ solution cl 	outp=prediccions; 
    

    random intercept TIME / subject=PATID type=UN g;
    
	ods output CovParms=params_cov; 
    
    title "Model 1: Linear Mixed-Effects Model";
run;

/* ============================================================== */

/* Multiple Imputation */

proc mi data=alzheimer25 out=mi_wide nimpute=10 seed=11 noprint;
    class SEX EDU JOB WZC;
    var SEX AGE EDU BMI JOB ADL WZC CDRSB0 bprs0-bprs6; 
    fcs nbiter=20; 
run;

proc sort data=mi_wide; by PATID _Imputation_; run;

data mi_long_final;
    set mi_wide;
    if _n_=1 then set stats; 

    array imp_arr[0:6] bprs0-bprs6;   

    do TIME = 0 to 6;
        BPRS = imp_arr[TIME];

        AGE_STD = (AGE - mean_age) / sd_age;
        BMI_STD = (BMI - mean_bmi) / sd_bmi;
        
        output;
    end;
    keep _Imputation_ PATID TIME BPRS SEX AGE AGE_STD EDU BMI BMI_STD JOB ADL WZC CDRSB0 TRIAL;
run;

proc sort data=mi_long_final;
    by _Imputation_ PATID;
run;

proc mixed data=mi_long_final method=REML plots(MAXPOINTS=None);
    by _Imputation_; 
    class PATID TRIAL SEX EDU JOB WZC / ref=first;
    

    model BPRS = TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
                 TIME*AGE TIME*EDU TIME*JOB TIME*ADL TIME*CDRSB0 / solution;
                 
    random intercept TIME / subject=PATID type=UN;
    ods output SolutionF=mixparms; 
    ods output CovParms=mi_cov_params;
run;

proc mianalyze parms=mixparms;
    class TRIAL SEX EDU JOB WZC;
    modeleffects Intercept TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
                 TIME*AGE TIME*EDU TIME*JOB TIME*ADL TIME*CDRSB0;
    title "LME Multiple Imputation";
run;

/* ============================================================== */

/* OLS residuals */
proc glm data=alzheimer_long_centered;
    class TRIAL SEX EDU JOB WZC;
    model BPRS = TIME TRIAL SEX AGE EDU BMI JOB ADL WZC CDRSB0
                 TIME*AGE TIME*EDU TIME*JOB TIME*ADL TIME*CDRSB0;
    output out=residus_ols r=resid_ols; 
run;

data _null_;
    set params_cov;
    name = compress(upcase(CovParm));
    if name = 'UN(1,1)'  then call symputx('d11', Estimate);
    if name = 'UN(2,1)'  then call symputx('d12', Estimate);
    if name = 'UN(2,2)'  then call symputx('d22', Estimate);
    if name = 'RESIDUAL' then call symputx('sigma2', Estimate);
run;

data funcio_teorica;
    do TIME = 0 to 6 by 0.1;
        fitted_marginal_var = &d11 + 2*TIME*&d12 + (TIME**2)*&d22 + &sigma2;
        output;
    end;
run;

data plot_final;
    set residus_ols funcio_teorica; 
    if resid_ols ne . then sq_resid_ols = resid_ols**2;
run;

proc sgplot data=plot_final;
    title "Comparison: LME Fitted Variance vs. OLS Empirical Variance";
    
    /* OLS Residuals */
    pbspline x=TIME y=sq_resid_ols / nomarkers lineattrs=(color=blue thickness=2) 
             legendlabel="Empirical Variance (OLS Residuals)";
    
    /* Fitted Variance */
    series x=TIME y=fitted_marginal_var / 
           lineattrs=(color=red thickness=3 pattern=dash) 
           legendlabel="Fitted Variance";
    
    xaxis label="Years" values=(0 to 6 by 1);
    yaxis label="Total Variance (BPRS)";
run;

/* ============================================================== */

proc means data=mi_cov_params noprint nway;
    class CovParm;
    var Estimate;
    output out=mi_avg_cov mean=avg_est;
run;

data _null_;
    set mi_avg_cov;
    name = compress(upcase(CovParm));
    if name = 'UN(1,1)'  then call symputx('d11_mi', avg_est);
    if name = 'UN(2,1)'  then call symputx('d12_mi', avg_est);
    if name = 'UN(2,2)'  then call symputx('d22_mi', avg_est);
    if name = 'RESIDUAL' then call symputx('sig2_mi', avg_est);
run;

data teorica_mi;
    do TIME = 0 to 6 by 0.1;
        fitted_var_mi = &d11_mi + 2*TIME*&d12_mi + (TIME**2)*&d22_mi + &sig2_mi;
        output;
    end;
run;

data plot_mi_final;
    set residus_ols teorica_mi;
    if resid_ols ne . then sq_resid_ols = resid_ols**2;
run;

proc sgplot data=plot_mi_final;
    title "MI Analysis: Averaged Fitted Variance vs. Empirical Variance";
    
    /* OLS Residuals */
    pbspline x=TIME y=sq_resid_ols / nomarkers lineattrs=(color=blue thickness=2) 
             legendlabel="Empirical Variance (OLS Residuals)";
    
    /* Fitted Variance MI */
    series x=TIME y=fitted_var_mi / lineattrs=(color=red thickness=3 pattern=dash) 
           legendlabel="Averaged Fitted Variance (MI)";
    
    xaxis label="Years" values=(0 to 6 by 1);
    yaxis label="Total Variance (BPRS)";
run;

/* ============================================================== */

data plot_total;
    set residus_ols    (keep=TIME resid_ols)
        funcio_teorica (keep=TIME fitted_marginal_var)
        teorica_mi     (keep=TIME fitted_var_mi);
    
    if resid_ols ne . then sq_resid_ols = resid_ols**2;
run;

proc sgplot data=plot_total;
    title "Variance Comparison: Empirical vs. LME vs. MI Models";
    
    /* OLS Residuals */
    pbspline x=TIME y=sq_resid_ols / nomarkers lineattrs=(color=blue thickness=2) 
             legendlabel="Empirical Variance (OLS Residuals)";
    
    /* Fitted Variance */
    series x=TIME y=fitted_marginal_var / lineattrs=(color=red thickness=2) 
           legendlabel="Fitted Variance (Original LME)";
    
    /* Fitted Variance MI */
    series x=TIME y=fitted_var_mi / lineattrs=(color=green thickness=3 pattern=dash) 
           legendlabel="Averaged Fitted Variance (MI Analysis)";
    
    xaxis label="Years" values=(0 to 6 by 1);
    yaxis label="Total Variance (BPRS)";
run;