*Data from Table 5.14: Space Shuttle disaster;
data shut1; input temp td @@; datalines;
 66 0 70 1 69 0 68 0 67 0 72 0 73 0 70 0 57 1 63 1
 70 1 78 0 67 0 53 1 67 0 75 0 70 0 81 0 76 0 79 0
 75 1 76 0 58 1
 ;
data shut2;
 do i=1 to 50; temp=i+29; td=.; output; end;
data shut3; set shut1 shut2;
proc logistic descending data=shut3; model td = temp / link=logit;
 output out=shut4 p=p1;
proc logistic descending data=shut3; model td = temp / link=cloglog;
 output out=shut5 p=p2;
proc logistic descending data=shut3; model td = temp / link=probit;
 output out=shut6 p=p3;
data shut7; set shut4 shut5 shut6;
proc sort data=shut7; by temp;
goptions;
 symbol1 color=black value=dot  interpol=none;
 symbol2 color=black value=none l=1 interpol=join;
 symbol3 color=black value=none l=2 interpol=join;
 symbol4 color=black value=none l=3 interpol=join;
 legend1 label=none value=('data' 'logit' 'cloglog' 'probit');
proc gplot data=shut7;
 plot td*temp p1*temp p2*temp p3*temp / overlay legend=legend1;
