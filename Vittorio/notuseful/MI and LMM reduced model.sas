data alzheimer25;
	set '/home/u64347574/alzheimer25.sas7bdat';
run;

proc means data=alzheimer25 noprint;
    var AGE BMI;
    output out=stats mean=mean_age mean_bmi std=sd_age sd_bmi;
run;

proc mi data=alzheimer25 seed=1234 out=alzheimer25_mi simple nimpute=10 round=1;
var trial edu inkomen job wzc age bmi sex adl taupet0 taupet1 taupet2 taupet3 taupet4 taupet5 taupet6
abpet0 abpet1 abpet2 abpet3 abpet4 abpet5 abpet6
cdrsb0 cdrsb1 cdrsb2 cdrsb3 cdrsb4 cdrsb5 cdrsb6
bprs0 bprs1 bprs2 bprs3 bprs4 bprs5 bprs6 ;
/*by patid;*/
run;

data alzheimer25_mi_base;
    set alzheimer25_mi;

    tau_base   = taupet0;
    ab_base    = abpet0;
    cdrsb_base = cdrsb0;
run;

data alzheimer_long;
    set alzheimer25_mi_base;
	
    array cdrsb_arr[0:6] cdrsb0-cdrsb6;
    array bprs_arr[0:6]  bprs0-bprs6;
    array abpet_arr[0:6] abpet0-abpet6;
    array taupet_arr[0:6] taupet0-taupet6;

    do TIME = 0 to 6;

        CDRSB  = cdrsb_arr[TIME];
        BPRS   = bprs_arr[TIME];
        ABPET  = abpet_arr[TIME];
        TAUPET = taupet_arr[TIME];

        /* CLASS time (se serve come categorica) */
        TIMECLSS = put(TIME, 1.);

        /* Una riga per soggetto-tempo */
        output;
    end;

    drop cdrsb0-cdrsb6
         bprs0-bprs6
         abpet0-abpet6
         taupet0-taupet6;
run;

proc means data=ALZHEIMER_LONG noprint;
    var ABPET TAUPET AGE BMI; 
    output out=bio_stats 
        mean=m_ab m_tau m_age m_bmi 
        std=s_ab s_tau s_age s_bmi;
run;

data alzheimer_long_centered;
    if _n_=1 then set bio_stats;
    set alzheimer_long;

    ABPET_STD = (ABPET - m_ab) / s_ab;
    TAUPET_STD = (TAUPET - m_tau) / s_tau;
    AGE_STD = (AGE - m_age) / s_age;
    BMI_STD = (BMI - m_bmi) / s_bmi;
run;

proc mixed data=alzheimer_long_centered method=REML;
    by _Imputation_;
    class sex patid;
    model bprs =  sex age bmi job
          adl wzc cdrsb_base
          tau_base
          time
          sex*time
          cdrsb_base*time
                  / solution;
    random intercept time / subject=patid type=un;
    ods output SolutionF=results;
run;

proc mianalyze parms=results;
    modeleffects Intercept sex age bmi job
          adl wzc cdrsb_base
          tau_base
          time
          sex*time
          cdrsb_base*time;
    ods output ParameterEstimates=mianalyze_results;
run;
