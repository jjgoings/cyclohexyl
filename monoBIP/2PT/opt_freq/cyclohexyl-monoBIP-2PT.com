%chk=cyclohexyl-monoBIP-2PT.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine monoBIP neutral, with both protons transferred
optimization

0 1
 C                 -0.52249400    2.72092200    0.33406500
 C                 -1.92037300    2.94681200    0.27176300
 C                 -2.46317700    4.22887600    0.21239200
 C                 -1.57039500    5.30477600    0.22296100
 C                 -0.18610600    5.10939800    0.28415400
 C                  0.37071300    3.81949400    0.33403700
 C                 -1.45275000    0.78025500    0.35284600
 H                 -3.53280100    4.38713500    0.16009000
 H                 -1.96040900    6.31387800    0.18207500
 H                  0.47727800    5.96642800    0.29018500
 N                 -2.47100700    1.68795100    0.28382800
 C                 -3.03159100   -1.19657800    0.29644300
 C                 -1.67158000   -0.65246700    0.39143600
 C                 -0.59527200   -1.51003600    0.52076200
 C                 -0.76828400   -2.90789100    0.57346600
 C                 -2.07158100   -3.44604300    0.47261300
 C                 -3.18620200   -2.64923300    0.33506900
 H                  0.40257800   -1.09377800    0.59201800
 H                 -2.19168900   -4.52361600    0.51264200
 C                  0.42057100   -3.80108100    0.75021400
 C                 -4.56946200   -3.20657400    0.23426100
 O                 -4.04374900   -0.45716500    0.18861900
 N                 -0.25164900    1.38069200    0.38343100
 C                  1.78915900    3.64287900    0.36612700
 H                  2.42742400    4.51992000    0.36827300
 N                  2.35761700    2.47525200    0.39037900
 C                  3.81495800    2.23978700    0.40237400
 C                  4.44185500    2.62112800   -0.94626400
 C                  4.49320500    2.97279600    1.56534700
 H                  3.92240600    1.16341000    0.55083000
 C                  5.94769000    2.33167800   -0.92898800
 H                  4.26967600    3.68929000   -1.12204200
 H                  3.94744100    2.07107400   -1.75098800
 C                  5.99683200    2.67409700    1.57442200
 H                  4.34164300    4.05163100    1.44719700
 H                  4.02543700    2.67642300    2.50739300
 C                  6.64425200    3.04412600    0.23518400
 H                  6.39102500    2.62850300   -1.88347600
 H                  6.09972300    1.24989100   -0.83273100
 H                  6.46973700    3.21878800    2.39610200
 H                  6.15055600    1.60592300    1.76974300
 H                  7.70846100    2.79267400    0.24793100
 H                  6.57367200    4.12915500    0.08952300
 H                 -3.44180800    1.38841600    0.24227800
 H                  1.72611800    1.66509600    0.38122600
 C                 -5.47246714   -2.61753227    1.33390255
 H                 -6.01074989   -1.77967677    0.94255888
 H                 -6.16478405   -3.36342095    1.66442845
 H                 -4.86927986   -2.30022998    2.15875670
 C                 -4.55450823   -4.74102559    0.36400961
 H                 -5.55439129   -5.11490814    0.29089525
 H                 -3.95900054   -5.15949759   -0.42030652
 H                 -4.13974290   -5.01481553    1.31159010
 C                 -5.21679804   -2.82685517   -1.11049536
 H                 -5.72404278   -1.89007672   -1.01025114
 H                 -4.45823658   -2.74324400   -1.86049115
 H                 -5.91788704   -3.58341433   -1.39508709
 C                  1.13363561   -3.48439407    2.07793696
 H                  1.62973304   -4.36274969    2.43469843
 H                  1.85266740   -2.70784642    1.92024409
 H                  0.41394714   -3.16255045    2.80137679
 C                  1.45054528   -3.57029916   -0.37116584
 H                  2.04882001   -2.71507308   -0.13546331
 H                  2.07912926   -4.43153910   -0.46087619
 H                  0.93931770   -3.40393696   -1.29629825
 C                  0.01572188   -5.28691191    0.75213411
 H                 -0.65505995   -5.47276102    1.56479138
 H                 -0.46835830   -5.52587625   -0.17169619
 H                  0.88929280   -5.89446109    0.86464124

radii=bondi
dis
rep
cav









