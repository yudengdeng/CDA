data drug;
input g r a c m count @@;
datalines;
0 1 1 1 1 405 0 1 1 1 0 268
0 1 1 0 1  13 0 1 1 0 0 218
0 1 0 1 1   1 0 1 0 1 0  17
0 1 0 0 1   1 0 1 0 0 0 117
1 1 1 1 1 453 1 1 1 1 0 228
1 1 1 0 1  28 1 1 1 0 0 201
1 1 0 1 1   1 1 1 0 1 0  17
1 1 0 0 1   1 1 1 0 0 0 133
0 0 1 1 1  23 0 0 1 1 0  23
0 0 1 0 1   2 0 0 1 0 0  19
0 0 0 1 1   0 0 0 0 1 0   1
0 0 0 0 1   0 0 0 0 0 0  12
1 0 1 1 1  30 1 0 1 1 0  19
1 0 1 0 1   1 1 0 1 0 0  18
1 0 0 1 1   1 1 0 0 1 0   8
1 0 0 0 1   0 1 0 0 0 0  17
;

run;

proc genmod;
 class g r a c m;
 model count=g|r|a g|r|c g|r|m g|a|c g|a|m g|c|m r|a|c r|a|m r|c|m a|c|m / link=log dist=poi type3;
 run;

proc genmod;
 class g r a c m;
 model count=g|r g|a g|c g|m r|a r|c r|m a|c a|m c|m / link=log dist=poi type3;
 run;

proc genmod;
 class g r a c m;
 model count=g|a g|m r|a r|m a|c a|m c|m / link=log dist=poi type3;
 run;
