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
 datalines;
 1 1 18  1 2 171 2 1 5 2 2 99
;
run;
proc freq  order=data; weight count;
 format Treatment tc. Outcome $oc.;
  tables Treatment*Outcome / norow nocol nopercent expected;
  exact chisq or;
  run;
proc genmod  order=data;  class Treatment Outcome;
 model count = Treatment Outcome Treatment*Outcome / dist=poi;
 run;
