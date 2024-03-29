*exercise 9.2;
data MBTI_data;
  input EI $ SN $ TF $ JP $ smoke drink n;
  datalines;
  e s t j 13 10 77 
  e s t p 11  8 42 
  e s f j 16  5 106 
  e s f p 19  7 79 
  e n t j 6 3 23 
  e n t p 4 2 18 
  e n f j 6 4 31 
  e n f p 23 15 80 
  i s t j 32 17 140 
  i s t p 9 3 52 
  i s f j 34  6 138 
  i s f p 29  4 106 
  i n t j 4 1 13 
  i n t p 9 5 35 
  i n f j 4 1 31 
  i n f p 22  6 79
  ;
run;
proc genmod order=data;
class EI SN TF JP;
model n = EI|SN|TF EI|SN|JP EI|TF|JP SN|TF|JP / dist=poi link=log type3;
run;
proc genmod order=data;
class EI SN TF JP;
model n = EI|SN|TF EI|SN|JP EI|TF|JP / dist=poi link=log type3;
run;
proc genmod order=data;
class EI SN TF JP;
model n = EI|SN|TF EI|SN|JP / dist=poi link=log type3;
run;
proc genmod order=data;
class EI SN TF JP;
model n = EI SN TF JP EI|SN EI|TF EI|JP SN|TF SN|JP TF|JP EI|SN|JP / dist=poi link=log type3;
run;
proc genmod order=data;
class EI SN TF JP;
model n = EI|SN EI|TF EI|JP SN|TF SN|JP TF|JP / dist=poi link=log type3;
run;
*Mutual independence;
proc genmod order=data;
class EI SN TF JP;
model n = EI SN TF JP / dist=poi link=log type3;
run;
*Homogeneous association;
proc genmod order=data;
class EI SN TF JP;
model n = EI|SN EI|TF EI|JP SN|TF SN|JP TF|JP / dist=poi link=log type3;
run;
*All three way factor;
proc genmod order=data;
class EI SN TF JP;
model
