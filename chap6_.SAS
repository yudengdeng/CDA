/* Using Logistic Models to Test Conditional Independence */
data cmh2; 
input center $ treat y n; treat=abs(treat-2);
datalines;
a 1 11 36
a 2 10 37
b 1 16 20
b 2 22 32
c 1 14 19
c 2 7 19
d 1 2 16
d 2 1 17
e 1 6 17
e 2 0 12
f 1 1 11
f 2 0 10
g 1 1 5
g 2 1 9
h 1 4 6
h 2 6 7
;
run;
*Conditional independence; 
proc logistic data=cmh2; class center; model y/n=center; run;
*Common odds ratio (additive model);
proc logistic data=cmh2; class center; model y/n=treat center; run;
*Test of goodness of fit for additive model;
proc logistic data=cmh2; class center; model y/n=treat center; run;
*Heterogenous odds ratio;
proc logistic data=cmh2; class center; model y/n= treat center treat*center; run;
run;
