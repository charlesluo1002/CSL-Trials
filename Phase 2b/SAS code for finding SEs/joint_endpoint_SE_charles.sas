*** Program to Compute SEs for Alternative CSL Designs; 
options notes source source2 errors=20; 
options mprint symbolgen; 
ods graphics off; 

*** Specify inputs; 
%let Ngroup = 250; 	* N/Group for base case scenario - Final csv file includes other Ns;
%let seed = 37;  
%let sige_1 = 0.67; %let sige_2 = 0.817; %let sige_3 = 1.15; * Defines within subject eGFR variance; 
%let Tint_1 = 0.083333; %let Tint_2 = 0.166667;  			 * Interval between visits, in yrs;
%let miss_1 = 0.05; %let miss_2 = 0.12; * Rate of loss to fup - 0.12 accounts for artificial censoring at med discontinuation;  
%let AcutePeriod = 0.25;									 * Allow for 3 month acute effect;

data one;
** Define combinations of input parameters for study design and variances components    **; 
** Modify as needed for the designs you want to consider								**; 
do isigc = 0 to 0; 					** Governs between pt slope variance;				
do ibgfr = 1 to 2;  				** Governs baseline eGFR;
do isigemult = 1 to 1;				** Governs within pt slope variance;
do iTint  = 1 to 2;					** Governs interval between eGFRs;
do iyrs  = 0 to 4;					** Governs follow-up period; 
do imiss = 0 to 1; 
if imiss = 0 then miss=&miss_1; if imiss = 1 then miss=&miss_2;
config + 1;							** Index for parameter configuration; 
sigc = 4 + (isigc)*1; 				** Assume Slope SD of 4; 
siga = 1;           				** Assume 1 ml/min/1.73m2/yr added acute slope SD over chronic slope SD;
if ibgfr = 1 then bgfr = 60; 		** Consider mean baseline eGFR of 40 or 60;
if ibgfr = 2 then bgfr = 40; 		** Consider mean baseline eGFR of 40 to 60;
if isigemult=1 then sigemult = &sige_1;
if isigemult=2 then sigemult = &sige_2;
if isigemult=3 then sigemult = &sige_3;
sige = sqrt(bgfr * sigemult);		** Assume residual variance proportoinal to baseline GFR;
Yrs  = 1 + iyrs/2; 					** Max Follow-up Yrs;
Accrual    = 1.5; 					** Period (yrs) of uniform patient accrual; 
AccrualAdj = Accrual/2;  			** Median time (yrs) of accrual;
YrChronic  = Yrs - &AcutePeriod - AccrualAdj;   ** Median Time Period for Chronic Slope;
YrTotal    = Yrs - AccrualAdj; 					** Median Time Period for Treatment; 
if iTint = 1 then Tint=&Tint_1; if iTint = 2 then Tint=&Tint_2; 
freq = int(.0001 + 1/Tint); ** Number eGFR measurements per year; 
*****************************************************************************************;
*****************************************************************************************;

**  Define design matrix based on input parameters                                     **;
**  The outcome y must be non-zero, but the actual value is ignored     			   **; 
do treatment = 0 to 1;
subno=0;
do isub = 1 to &Ngroup;
  id+1;subno+1;
  isubYrs   = yrs - Accrual*sqrt(isub/&Ngroup);                  ** Follow-up period for subject isub;
  isubChron = yrs - Accrual*sqrt(isub/&Ngroup) - &AcutePeriod;   ** Duration of chronic phase for subject isub; 
  t = 0; y = rannor(&seed);output; 				** add one extra baseline eGFR;
  do it = 0 to int(freq*isubyrs); 
    t = it*Tint; 								** follow-up time; 
    y = rannor(&seed); 							** random value for outcome, will be ignored; 
	output; 
	end; 		
 t = &AcutePeriod; y = rannor(&seed); 		    ** add 2nd eGFR at end of acute phase; 
 output; 
 t = isubyrs; y = rannor(&seed);    			** end closeout eGFR;
 output;  
 t = isubyrs+&Acuteperiod; y = rannor(&seed);output;   ** add off-treatment eGFR; 
 t = isubyrs+&Acuteperiod; y = rannor(&seed);output;   ** add second off-treatment eGFR;
 end; end;                     					** end subject and treatment group loops; 
 end; end; end; end; end; end;         			** end parameter configuration loops; 
run; 

data oneb; set one; 
** Define X and Z matrices for linear mixed model;
** X is fixed effect design matrix;
** Z is random effect design matrix; 
int = 1; 
timeA = min(0,t-&AcutePeriod);  			** 1st term for 3-slope linear spline in time; 
timeB = max(0,t-&AcutePeriod);				** 2nd term for 3-slope linear spline in time; 
timeB = min(timeB,isubChron);				** 2nd term for 3-slope lienar spline in time; 
timeC = max(0,t-(isubChron+&Acuteperiod));	** 3rd term for 3-slope linear spline in time; 
timeBRan = t-&AcutePeriod;

******************************* Insert missing data ***********************************; 
** First generate alternate sequence of subject nubmers that is approximately orthogonal to original subno;
submod2 =  subno - 2*int(subno/2); 
submod4 =  subno - 4*int(subno/4); 
submod8 =  subno - 8*int(subno/8); 
submod16 = subno - 16*int(subno/16);
submod32 = subno - 32*int(subno/32); 
proc sort data = oneb; by isigc ibgfr isigemult iTint iyrs imiss treatment submod2 submod4 submod8 submod16 submod32 subno t;
data onec; set oneb;   by isigc ibgfr isigemult iTint iyrs imiss treatment submod2 submod4 submod8 submod16 submod32 subno t; 
if first.treatment then altsubno = 0;
if first.subno then altsubno+1;
retain; 
proc sort data=onec;   by isigc ibgfr isigemult iTint iyrs imiss treatment subno t;
run;
data onec; set onec; 
ptmiss = 1 - exp(-t*&miss);
ntmiss = int(&Ngroup*ptmiss);
if altsubno <= ntmiss then y=.;
run;

%macro mixed;
data psumall; set _NULL_;
*** the lines below that define the looping variables much exactly match the same loops at the beginning of the program above;
%do isigc = 0 %to 0;
%do ibgfr = 1 %to 2;
%do isigemult = 1 %to 1;
%do iTint  = 1 %to 2;
%do iyrs  = 0 %to 4;
%do imiss = 0 %to 1;

data indat; set onec; 
 if isigc = &isigc;
 if ibgfr = &ibgfr;
 if isigemult = &isigemult;
 if iTint  = &iTint;
 if iyrs  = &iyrs;
 if imiss = &imiss;
data first; set indat; by config; if first.config; run;

data _NULL_; set first; 
call symput ('SIGA',put(siga,7.4));
call symput ('SIGC',put(sigc,7.4));
call symput ('SIGE',put(sige,7.4));
call symput ('SIGEMULT',put(sigemult,7.4));
call symput ('BGFR',put(bgfr,7.4));
call symput ('Tint',put(Tint,7.4));
call symput ('YrTotal',put(yrs,7.4));
call symput ('YrChronic',put(yrchronic,7.4)); 
call symput ('MissYr',put(miss,7.4));
run;

proc contents data=first; run;

%macro VarComp;
%global vara varc vare;
%let vara = %sysevalf(&SIGA * &SIGA); 
%let varc = %sysevalf(&SIGC * &SIGC);
%let vare = %sysevalf(&SIGE * &SIGE);
%mend;
run;
%VarComp;
run;

** run proc mixed with noiter option to obtain SEs for treatment effects for Chronic, Total and Offtreatment slopes;
ods select estimates; ods output estimates=estout; 
proc mixed data=indat noprofile;
 model y = timeA timeB timeC treatment*timeA treatment*timeB treatment*timeC;
 random intercept timeA timeBRan timeC/sub=id type=un;
 parms (100) (0) (&VARA) (0) (0) (&VARC) (0) (0) (0) (&VARA) (&VARE)/ noiter;
 *** obtain SEs for estimated tretment effects on relevant slope parameters;  
 estimate 'Chronic' treatment*timeB 1;								
 estimate 'Total'   treatment*timeA &AcutePeriod treatment*timeB &YrChronic/divisor=&YrTotal;
 estimate 'OffTrt'  treatment*timeA &AcutePeriod treatment*timeB &YrChronic treatment*timeC &AcutePeriod/divisor=&YrTotal; 
 run;

data psum; set estout;
siga = &SIGA; sigc=&SIGC; sige=&SIGE; sigemult=&SIGEMULT; bgfr = &BGFR; YrTotal = &YrTotal; Tint=&Tint;MissYr=&MissYr; 
YrChronic = &YrChronic;
 keep siga sigc sige sigemult bgfr Tint YrTotal YrChronic MissYr label stderr; 
 run;
data psumall; set psumall psum;
proc datasets; delete estout psum; run; quit; 
%end; %end; %end; %end; %end; %end;
%mend; 
run;
%mixed; 
run;
options notes source source2 errors=20; 
options mprint symbolgen;  run;


proc freq data = psumall; table YrFup; run;

*** estimate SEs for different sample sizes;
data n150; 
set psumall; ngroup = 150; stderr = stderr*sqrt(250/150);
data n200; 
set psumall; ngroup = 200; stderr = stderr*sqrt(250/200);
data n250; 
set psumall; ngroup = 250; stderr = stderr*sqrt(250/250);
data n1000; 
set psumall; ngroup = 1000; stderr = stderr*sqrt(250/1000);
data ninf; 
set psumall; ngroup = 100000000; stderr = stderr*sqrt(250/100000000);
data estsum; set n150 n200 n250 n1000 ninf;  
if label = "Chronic";    *** This line is optional: Use only if you want to restrict attention to the chronic slope;
proc export data=estsum
   outfile='H:\My Documents\CSL work file\further work\completeSlopeSEs.csv'
   dbms=csv
   replace;
run;
