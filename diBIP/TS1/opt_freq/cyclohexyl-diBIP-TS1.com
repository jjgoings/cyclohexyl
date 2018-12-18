%chk=cyclohexyl-diBIP-TS1.chk
%mem=10GB
%nprocshared=16
#p opt=(readfc,recalcfc=6,tight,ts,noeigentest) freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm guess=read
empiricaldispersion=gd3bj int=grid=ultrafine

cyclohexylimine diBIP oxidized, first proton transition state optimization

1 2
 C                  1.67485500    1.39685500    0.15875300
 C                  0.45604500    2.09591300   -0.00246200
 C                  0.45803200    3.49823300   -0.00670400
 C                  1.67772900    4.14919500    0.16296600
 C                  2.87813100    3.43524800    0.33384300
 C                  2.91164700    2.03663500    0.33349600
 C                 -0.01130700   -0.01514700   -0.06204000
 H                 -0.46139200    4.05468300   -0.14002400
 H                  1.70758300    5.23242000    0.16573700
 H                  3.80775100    3.97754900    0.46896400
 C                 -2.15435000   -1.24570900   -0.13087100
 C                 -2.94038700   -2.42942900   -0.16496800
 C                 -2.32064500   -3.68908900   -0.20948000
 C                 -0.93441300   -3.71158000   -0.22230900
 C                 -0.16187900   -2.52894300   -0.17993600
 C                 -0.74786900   -1.26475400   -0.13340800
 C                 -4.29849000   -0.73071200   -0.08183600
 H                 -2.90583300   -4.59886400   -0.23533900
 H                 -0.41863700   -4.66281100   -0.26342300
 H                  0.91826300   -2.61241700   -0.19032200
 N                 -0.58305500    1.18579700   -0.13392800
 C                 -6.76867000   -0.76373400   -0.23318500
 C                 -5.53864200   -0.00463300   -0.06326900
 C                 -5.61964800    1.37534300    0.08435100
 C                 -6.85513100    2.02716400    0.04763000
 C                 -8.04900900    1.26458600   -0.11968100
 C                 -8.04105300   -0.09901100   -0.25859500
 H                 -4.71948700    1.95562600    0.23762700
 H                 -8.99446900    1.79640400   -0.12791600
 C                 -6.96408300    3.50727400    0.19061700
 H                 -7.59762600    3.74336800    1.05240300
 H                 -5.99626100    3.99158000    0.32736300
 H                 -7.46751500    3.93751900   -0.68075600
 C                 -9.28468200   -0.91895300   -0.39854600
 H                 -9.16784100   -1.68753800   -1.16484300
 H                 -9.50187700   -1.43398100    0.54291100
 H                -10.13457300   -0.28500100   -0.65015200
 O                 -6.75550800   -2.05955700   -0.35745300
 N                 -3.03587400   -0.20204200   -0.07556000
 N                 -4.26325500   -2.07232000   -0.12650300
 H                 -2.72007700    0.76568900   -0.09727000
 N                  1.34698100    0.06551900    0.10974800
 H                  2.03925200   -0.66568200    0.23151700
 C                  4.15422000    1.28724000    0.51140700
 H                  5.05636100    1.89281100    0.62743900
 N                  4.16316600    0.00689000    0.52912500
 C                  5.42433800   -0.72304000    0.66940300
 C                  5.62203100   -1.57707000   -0.59745900
 C                  6.68523000    0.09429000    0.97057500
 H                  5.27168900   -1.42164000    1.50421200
 C                  6.84436200   -2.49128500   -0.47020600
 H                  5.74986000   -0.89685800   -1.44880900
 H                  4.71348600   -2.15714800   -0.78280600
 C                  7.89578600   -0.83585900    1.12171800
 H                  6.87049400    0.80026600    0.15130100
 H                  6.53909900    0.68604500    1.88049400
 C                  8.10481700   -1.68523800   -0.13711000
 H                  6.98835500   -3.05604900   -1.39689000
 H                  6.66608500   -3.22285100    0.32751000
 H                  8.79427100   -0.25464900    1.34752300
 H                  7.73257500   -1.50604700    1.97522200
 H                  8.95870200   -2.35609200   -0.00258200
 H                  8.34672000   -1.02713800   -0.98122600
 H                 -5.76540500   -2.38491900   -0.29965100


radii=bondi
dis
rep
cav





