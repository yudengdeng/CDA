proc format;
value lakefmt 1='Hancock'
              2='Oklawaha'
              3='Trafford'
              4='George';
value gdrfmt  1='Male'
              2='Female';
value foodfmt 1='Fish'
              2='Invertebrate'
              3='Reptile'
              4='Bird'
              5='Other';
value sizefmt 1='<=2.3 m'
              2='>2.3 m';
data alligator;
input Lake Gender Size Food Count ;
label Lake="Lake" Gender="Gender" Size="Size (m)" Food="Primary Food Choice";
format Lake lakefmt. Gender gdrfmt. Food foodfmt. Size sizefmt.;
cards;
1 1 1 1 7   
1 1 1 2 1   
1 1 1 3 0   
1 1 1 4 0   
1 1 1 5 5
1 1 2 1 4   
1 1 2 2 0   
1 1 2 3 0   
1 1 2 4 1   
1 1 2 5 2        
1 2 1 1 16  
1 2 1 2 3   
1 2 1 3 2   
1 2 1 4 2   
1 2 1 5 3
1 2 2 1 3   
1 2 2 2 0   
1 2 2 3 1   
1 2 2 4 2   
1 2 2 5 3
2 1 1 1 2   
2 1 1 2 2   
2 1 1 3 0   
2 1 1 4 0   
2 1 1 5 1  
2 1 2 1 13  
2 1 2 2 7   
2 1 2 3 6   
2 1 2 4 0   
2 1 2 5 0  
2 2 1 1 3   
2 2 1 2 9   
2 2 1 3 1   
2 2 1 4 0   
2 2 1 5 2  
2 2 2 1 0   
2 2 2 2 1   
2 2 2 3 0   
2 2 2 4 1   
2 2 2 5 0 
3 1 1 1 3   
3 1 1 2 7   
3 1 1 3 1   
3 1 1 4 0   
3 1 1 5 1
3 1 2 1 8   
3 1 2 2 6   
3 1 2 3 6   
3 1 2 4 3   
3 1 2 5 5  
3 2 1 1 2   
3 2 1 2 4   
3 2 1 3 1   
3 2 1 4 1   
3 2 1 5 4          
3 2 2 1 0   
3 2 2 2 1   
3 2 2 3 0   
3 2 2 4 0   
3 2 2 5 0  
4 1 1 1 13  
4 1 1 2 10  
4 1 1 3 0   
4 1 1 4 2   
4 1 1 5 2  
4 1 2 1 9   
4 1 2 2 0   
4 1 2 3 0   
4 1 2 4 1   
4 1 2 5 2  
4 2 1 1 3   
4 2 1 2 9   
4 2 1 3 1   
4 2 1 4 0   
4 2 1 5 1
4 2 2 1 8   
4 2 2 2 1   
4 2 2 3 0   
4 2 2 4 0   
4 2 2 5 1    
;
run;       
/* Reproduce format of Agresti's table   */
proc tabulate data=alligator format=f12.0;
keylabel sum=' ';
class lake gender size food; var count;
table lake*gender*size, food*count=' '/rts=30; *Row title space;
run;
*Tim's backward elimination step;
/* The final category is the default baseline--use of REF changes that */
proc logistic data=alligator; 
freq count;
class lake(ref='George') size gender / param=ref;
/* You can use PARAM=EFFECT to obtain catmod-like coding */
model food(ref='Fish') = lake size gender lake*size size*gender  lake*gender / link=glogit 
  aggregate=(lake size gender) scale=none selection=backward;
/* AGGREGATE and SCALE are needed in tandem to obtain fit statistics */
/* AGGREGATE works much as POPULATION does in catmod */
run;
*The reduced model;
proc logistic data=alligator; 
freq count;
class lake(ref='George') size gender / param=ref;
/* You can use PARAM=EFFECT to obtain catmod-like coding */
model food(ref='Fish') = lake size / link=glogit 
  aggregate=(lake size) scale=none;
run;
*Similar model in PROC CATMOD, but with effects-style coding;
proc catmod data=alligator;
weight count;
/* We need the POPULATION statement so that reduced models are hierarchical */
population lake size;
model food = lake size / noiter;
run;

*Reproducing the CATMOD output in GENMOD;
proc logistic data=alligator; 
freq count;
class lake(ref='George') size gender / param=effect;
/* You can use PARAM=EFFECT to obtain catmod-like coding */
model food(ref='Other') = lake size / link=glogit 
  aggregate=(lake size gender) scale=none;
run;

*Birds as reference category;
proc logistic data=alligator; 
freq count;
class lake(ref='George') size gender / param=ref;
model food(ref='Bird') = lake size / link=glogit 
  aggregate=(lake size) scale=none;
run;

*Reproducing the CATMOD output in GENMOD;
proc logistic data=alligator; 
freq count;
class lake(ref='George') size gender / param=effect;
/* You can use PARAM=EFFECT to obtain catmod-like coding */
model food(ref='Other') = lake size / link=glogit 
  aggregate=(lake size gender) scale=none;
run;
