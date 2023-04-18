data ideology;
*Democratic party = 1, Republican party = 0
*Very Liberal ideology = 1, Slightly Liberal ideology = 2, Moderate ideology = 3, Slightly Conservative ideology = 4, Very Conservative ideology = 5
;
input party ideology count @@;
cards;
1 1 80   
1 2 81   
1 3 171   
1 4 41   
1 5 55   
0 1 30   
0 2 46   
0 3 148   
0 4 84   
0 5 99
;
proc catmod data = ideology;
  weight count; response alogits;
  model ideology = party;
run;
proc catmod data = ideology;
 weight count; response alogits;
 model ideology = _response_ party;
run;
