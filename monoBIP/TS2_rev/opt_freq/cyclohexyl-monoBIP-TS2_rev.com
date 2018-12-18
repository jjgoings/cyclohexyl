%chk=cyclohexyl-monoBIP-TS2_rev.chk
%mem=10GB
%nprocshared=16
#p opt=(calcfc,tight,ts,noeigentest) freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm guess=read
empiricaldispersion=gd3bj int=grid=ultrafine

cyclohexylimine monoBIP oxidized, optimize to transition state of first
proton transferred

1 2
 C                 -0.42419100    2.67174000    0.27773200
 C                 -1.83568200    2.79278600    0.23762000
 C                 -2.43393200    4.05878100    0.14154700
 C                 -1.58508000    5.15590000    0.09356900
 C                 -0.17867900    5.02067400    0.13607900
 C                  0.44283200    3.77616900    0.22938600
 C                 -1.33668500    0.66603600    0.36926400
 H                 -3.51006600    4.16886900    0.10552300
 H                 -2.00792600    6.15043800    0.02130800
 H                  0.43820000    5.91088900    0.09588900
 N                 -2.36559100    1.52550800    0.29775000
 C                 -2.95818400   -1.18914800    0.41953000
 C                 -1.56604400   -0.75274500    0.42968000
 C                 -0.55905600   -1.70280700    0.47879500
 C                 -0.86680700   -3.07204700    0.52101300
 C                 -2.23064000   -3.48478300    0.52576500
 C                 -3.27348400   -2.59449300    0.48095200
 H                  0.47865100   -1.39238600    0.47916500
 H                 -2.44452100   -4.54683600    0.56929100
 C                  0.21415300   -4.10109200    0.53956400
 H                  1.19640100   -3.65187800    0.68097600
 H                  0.03505000   -4.83794700    1.32772300
 H                  0.21665400   -4.65346800   -0.40708400
 C                 -4.70983400   -3.00981100    0.48526000
 H                 -5.24634700   -2.54246200    1.31532000
 H                 -5.20862100   -2.69063800   -0.43430000
 H                 -4.79374200   -4.09257600    0.57394700
 O                 -3.92255000   -0.33374200    0.35053100
 N                 -0.04176810    1.37848795    0.36166845
 H                  1.31235115    1.59084831    0.38055655
 C                  1.89445000    3.60335200    0.27084400
 H                  2.49113600    4.51934300    0.21673600
 N                  2.30150610    2.38779005    0.36581755
 C                  3.76569210    2.26830705    0.40533555
 C                  4.41443310    2.69491705   -0.92192645
 C                  4.39559310    3.03243005    1.58144755
 H                  3.97548810    1.20273005    0.55061155
 C                  5.93004010    2.47003805   -0.87863245
 H                  4.20120710    3.75745805   -1.08922145
 H                  3.96128810    2.13861605   -1.74792145
 C                  5.90999810    2.80109705    1.62865855
 H                  4.19807810    4.10389605    1.45798855
 H                  3.91964010    2.72242105    2.51632155
 C                  6.57071110    3.20382705    0.30488755
 H                  6.38471710    2.79095805   -1.82061745
 H                  6.12658810    1.39453805   -0.78435145
 H                  6.34520010    3.36198505    2.46137355
 H                  6.10675710    1.73945105    1.82366055
 H                  7.64583810    3.00346105    0.33957155
 H                  6.45294410    4.28535305    0.16180755
 H                 -3.45550900    0.69606700    0.30997700

radii=bondi
dis
rep
cav











