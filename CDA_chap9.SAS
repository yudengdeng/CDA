*loglinear analysis of I by J table;
proc format;
value ic 1='Twisted' 2='Direct blow' 3='Both';
value rc 1='Excellent' 2='Good' 3='Fair-to-poor';
run;
data table;
 input Injury Result count @@;
 format Injury ic. Result rc.;
 datalines;
1 1 21 1 2 11 1 3 4 2 1 3 2 2 2 2 3 2 3 1 7 3 2 1 3 3 1
;
run;
proc freq data=table order=data; weight count;
 tables Injury*Result / chisq;
 run;
proc genmod data=table order=data;  class Injury Result;
 model count = Injury Result / dist=poi;
 run;
*loglinear analysis of 2 by 2 table;
proc format;
  value tc 1='Placebo' 2='Aspirin';
  value oc 1='Fatal' 2='Nonfatal';
run;
data table;
  input Treatment Outcome count @@;
  format Treatment $tc. Outcome $oc.;
  datalines;
 1 1 18  1 2 171 2 1 5 2 2 99
;
run;
proc freq data=table order=data; weight count;
  tables Treatment*Outcome / norow nocol nopercent expected;
  exact chisq or;
run;
proc genmod  order=data;  class Treatment Outcome;
*Seat belt;
data belt;
input g l b i count @@; 
datalines;
1 1 1 1  7287  
1 1 1 2   996   
1 1 2 1 11587  
1 1 2 2   759
1 2 1 1  3246  
1 2 1 2   973   
1 2 2 1  6134  
1 2 2 2   757
2 1 1 1 10381  
2 1 1 2   812   
2 1 2 1 10969  
2 1 2 2   380
2 2 1 1  6123  
2 2 1 2  1084   
2 2 2 1  6693  
2 2 2 2   513
;
run;

proc genmod order=data; class g l b i;
 model count = g|l|b g|l|i l|b|i g|b|i / dist=poi link=log type3;
run;

proc genmod order=data; class g l b i;
 model count = g|l|b g|l|i l|b|i  / dist=poi link=log type3;
run;


proc genmod order=data; class g l b i;
 model count = g|l|b l|b|i g|i / dist=poi link=log type3;
run;

******************************************************
Collapsing over B and L affects the relationship
between G and I??
******************************************************;
proc genmod order=data; class g l b i;
 model count = g|i  / dist=poi link=log type3;
run;
 model count = Treatment Outcome Treatment*Outcome / dist=poi;
run;

data drugs;
input pers chol bp count @@;
datalines;
1 1 1 716 1 1 2 79
1 2 1 207 1 2 2 25
2 1 1 819 2 1 2 67
2 2 1 186 2 2 2 22
;
proc genmod order=data; class pers chol bp;
 model count = pers|chol|bp / dist=poi link=log type3;
run;
proc genmod order=data; class pers chol bp;
 model count = pers|chol pers|bp chol|bp / dist=poi link=log type3;
run;
