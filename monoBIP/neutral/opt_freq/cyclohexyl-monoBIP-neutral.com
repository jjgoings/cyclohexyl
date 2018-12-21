%chk=cyclohexyl-monoBIP-neutral.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine monoBIP neutral optimization

0 1
 C                 -0.46401600    2.68790700    0.30291800
 C                 -1.86559200    2.86169400    0.26695800
 C                 -2.40612800    4.15085300    0.19417700
 C                 -1.51927000    5.22726300    0.16518300
 C                 -0.12823100    5.03780000    0.20066400
 C                  0.44092000    3.75929400    0.26813400
 C                 -1.48608500    0.71971000    0.36879800
 H                 -3.47818600    4.30209700    0.16241700
 H                 -1.91162100    6.23610100    0.11298800
 H                  0.53032600    5.89928800    0.17620900
 N                 -2.47353400    1.61652300    0.30846400
 C                 -3.02750500   -1.21427700    0.36500200
 C                 -1.70217200   -0.71610200    0.41777600
 C                 -0.62863100   -1.61999900    0.49939800
 C                 -0.83819100   -2.99414800    0.53310200
 C                 -2.16213200   -3.45784500    0.47675700
 C                 -3.25718600   -2.60175000    0.39629500
 H                  0.38567300   -1.23676900    0.54236300
 H                 -2.34228500   -4.52854100    0.50125800
 C                  0.31600300   -3.96056300    0.62940000
 C                 -4.66791900   -3.11878400    0.33948400
 O                 -4.10482100   -0.40263400    0.28475800
 H                 -3.77383200    0.54022700    0.28035300
 N                 -0.25751000    1.33090900    0.36813600
 H                  0.67240900    0.92968100    0.39508100
 C                  1.88554300    3.55256300    0.29463800
 H                  2.49516200    4.46137000    0.24061200
 N                  2.39390200    2.37861800    0.37762400
 C                  3.85928400    2.26296700    0.40382400
 C                  4.49850600    2.69909400   -0.92481900
 C                  4.49839900    3.01984000    1.58009400
 H                  4.07513600    1.19715500    0.54110600
 C                  6.01567500    2.48271100   -0.89580100
 H                  4.27937400    3.76153000   -1.08518800
 H                  4.04171800    2.14564900   -1.75075800
 C                  6.01426100    2.79581500    1.61345600
 H                  4.29537900    4.09163800    1.46792200
 H                  4.03104900    2.70108500    2.51648700
 C                  6.66129000    3.21274300    0.28743300
 H                  6.46151900    2.81246700   -1.83920300
 H                  6.21981600    1.40781900   -0.80990300
 H                  6.45458800    3.35172300    2.44704700
 H                  6.21820500    1.73349300    1.79721900
 H                  7.73840300    3.02073200    0.31144500
 H                  6.53433900    4.29440000    0.15274700
 C                  1.35109905   -3.71903866   -0.48497972
 H                  2.04715567   -2.97030709   -0.16905853
 H                  1.87478657   -4.62959519   -0.68878893
 H                  0.85054467   -3.38940141   -1.37136824
 C                  1.05667327   -3.82032991    1.97228579
 H                  0.34406826   -3.71046874    2.76287078
 H                  1.64966790   -4.69397997    2.14547053
 H                  1.69090500   -2.95910636    1.94155967
 C                 -0.17789242   -5.41470788    0.51481379
 H                 -0.84186669   -5.63080386    1.32557929
 H                 -0.69440921   -5.54386265   -0.41331845
 H                  0.65943774   -6.07980413    0.55256548
 C                 -5.35968348   -2.71715829   -0.97647903
 H                 -5.97516069   -3.52368768   -1.31647630
 H                 -5.96591743   -1.85115626   -0.81088652
 H                 -4.61861362   -2.49757941   -1.71641185
 C                 -5.52062827   -2.53414666    1.48083674
 H                 -6.02247195   -1.65414400    1.13638320
 H                 -6.24401135   -3.25832327    1.79258208
 H                 -4.88786806   -2.28376382    2.30656245
 C                 -4.68987798   -4.65458692    0.45096130
 H                 -4.15043671   -5.07950319   -0.36961734
 H                 -4.23188438   -4.95198590    1.37112200
 H                 -5.70257008   -4.99935508    0.42883420

radii=bondi
dis
rep
cav





