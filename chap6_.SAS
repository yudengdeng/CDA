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
/* CMH Testing*/
data cmh;
 input center $ treat response count;
 datalines;
 a 1 1 11
 a 1 2 25
 a 2 1 10
 a 2 2 27
 b 1 1 16
 b 1 2 4
 b 2 1 22
 b 2 2 10
 c 1 1 14
 c 1 2 5
 c 2 1 7 
 c 2 2 12
 d 1 1 2
 d 1 2 14
 d 2 1 1
 d 2 2 16
 e 1 1 6 
 e 1 2 11
 e 2 1 0 
 e 2 2 12
 f 1 1 1
 f 1 2 10
 f 2 1 0 
 f 2 2 10
 g 1 1 1
 g 1 2 4
 g 2 1 1
 g 2 2 8
 h 1 1 4
 h 1 2 2
 h 2 1 6
 h 2 2 1
;
run;
proc freq data=cmh; 
weight count; 
tables center*treat*response / cmh;
run;
