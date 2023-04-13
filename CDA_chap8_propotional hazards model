data leuk1;
 input x1 x2 y @@;
 datalines;
1    6.62 3 1    7.74 2 1    8.36 2 1    7.86 3 1    8.69 1 1    9.25 3
1    9.21 3 1    9.74 1 1    8.59 1 1    8.85 3 1    9.14 2 1   10.37 1
1   10.46 1 1   10.85 1 1   11.51 1 1   11.51 1 1   11.51 2 0    8.38 2
0    8.00 2 0    8.29 1 0    7.31 1 0    9.10 1 0    8.57 1 0    9.21 1
0    9.85 1 0   10.20 1 0   10.23 1 0   10.34 1 0   10.16 1 0    9.95 1
0   11.27 1 0   11.51 1 0   11.51 1
;
run;
proc nlmixed; * effect of beta constant across stages;
 parms a1=-7 a2=-6 b1=-3 b2=1; * started with a1=0 a2=1 b1=0 b2=0;
 p1=exp(a1+x1*b1+x2*b2); p2=exp(a2+x1*b1+x2*b2);
 if (y=1) then z=(p1/(1+p1));
 if (y=2) then z=(1/(1+p1))*(p2/(1+p2));
 if (y=3) then z=(1/(1+p1))*(1/(1+p2));
 if (z>1e-8) then ll=log(z); else ll=-1e100;
 model y ~ general(ll);
 run;
 *Proportional odds;
 proc logistic data=leuk1; 
 model y=x1 x2;
 run;
