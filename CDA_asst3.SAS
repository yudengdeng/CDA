*Table 5.14 from Agresti 3rd;
data shuttle_data;
   input Ft Temp TD;
   datalines;
   1 66 0
   6 72 0
   11 70 1
   16 75 0
   21 75 1
   2 70 1
   7 73 0
   12 78 0
   17 70 0
   22 76 0
   3 69 0
   8 70 0
   13 67 0
   18 81 0
   23 58 1
   4 68 0
   9 57 1
   14 53 1
   19 76 0
   5 67 0
   10 63 1
   15 67 0
   20 79 0
   ;
run;
proc logistic data=shuttle_data outest=shuttle_data_coefficients;
	class TD(ref='0');
	model TD = Temp;
	output out=preds predicted=PredProb;
run;
proc print data = shuttle_data_coefficients;
run;
/* Load the saved coefficient estimates */
data prediction;
	set shuttle_data_coefficients;
	logit = Intercept + Temp*31;
	prob = exp(logit) / (1 + exp(logit));
	put prob; *read the output from Log;
run;
*figure of fitted lines;
proc sgplot data=preds;
   scatter x=Temp y=TD / markerattrs=(symbol=circlefilled);
   loess x=Temp y=PredProb / lineattrs=(color=blue);
   yaxis label="TD" values=(0 to 1 by 0.2);
   xaxis label="Temperature (Fahrenheit)";
run;
*Plot a variety of diagnostic plots;
proc logistic data=shuttle_data descending plots=all;
class TD(ref='0');
model TD = Temp;
output out=diag2 stdreschi=r c=c;
run;
proc sgscatter data=diag2;
title "Cook's distance and Std. Pearson resids";
plot c*Ft r*Ft;
run;
*Print unusual cases;
proc print data=diag2(where=(c>0.2 or r>3 or r<-3));
var TD Temp c r;
run;
*Table 10.1 from Agresti 3rd;
data Table_10_1;
input a c m r g count;
datalines;
1 1 1 1 1 405      
1 1 1 2 1  23      
1 2 1 1 1  13      
1 2 1 2 1   2      
2 1 1 1 1   1      
2 1 1 2 1   0      
2 2 1 1 1   1      
2 2 1 2 1   0      
1 1 2 1 1 268
1 1 2 2 1  23
1 2 2 1 1 218
1 2 2 2 1  19
2 1 2 1 1  17
2 1 2 2 1   1
2 2 2 1 1 117
2 2 2 2 1  12
1 1 1 1 2 453
1 1 1 2 2  30
1 2 1 1 2  28
1 2 1 2 2   1
2 1 1 1 2   1
2 1 1 2 2   1
2 2 1 1 2   1
2 2 1 2 2   0
1 1 2 1 2 228
1 1 2 2 2  19
1 2 2 1 2 201
1 2 2 2 2  18
2 1 2 1 2  17
2 1 2 2 2   8
2 2 2 1 2 133
2 2 2 2 2  17
;
run;
*Backward selection with slstay=0.05;
proc logistic data=Table_10_1 descending;
	class a c m r g/param=ref;
	WEIGHT count;
	model m = a c r g a*c a*r a*g c*r c*g r*g/selection=backward;
run;
*asst 3 problem 6;
data berkeley;
input dept$ gender$ admit not_admit @@;
total=admit+not_admit;
datalines;
a male 512 313 a female 89 19
b male 353 207 b female 17 8
c male 120 205 c female 202 391
d male 138 279 d female 131 244
e male 53 138 e female 94 299
f male 22 351 f female 24 317
;
*Additive model;
proc logistic data=berkeley; class dept gender / param=ref;
model admit/total=dept gender / aggregate scale=none;
*Model with intercation;
proc logistic data=berkeley; class dept gender / param=ref;
model admit/total=dept gender dept*gender / aggregate scale=none;
*LRT;
proc iml;
x = 877.056 - 856.85;    /* value to evaluate CDF at */
df = 11 - 6;    /* degrees of freedom */
prob = 1 - cdf("chisquare", x, df);    /* calculate CDF */
print prob;
quit;
