%chk=cyclohexyl-triBIP-E0PT.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine triBIP oxidized, no proton transfer optimization

1 2
 C                  1.67909300    0.39251600    0.95170400
 C                  0.58313800    1.23902000    0.66246000
 C                  0.78508400    2.61740800    0.52156300
 C                  2.08217500    3.10051200    0.67602900
 C                  3.16190400    2.24961300    0.96754400
 C                  2.98961900    0.86785300    1.11401700
 C                 -0.19030100   -0.77716200    0.79391400
 H                 -0.04168800    3.28013900    0.29793300
 H                  2.27041100    4.16211200    0.56809900
 H                  4.15356500    2.67443900    1.06984500
 C                 -2.47690300   -1.67209000    0.54491100
 C                 -3.44440300   -2.71240100    0.53535400
 C                 -3.05805400   -4.03895800    0.78370800
 C                 -1.71293300   -4.26881900    1.03011200
 C                 -0.75691600   -3.22733900    1.03565700
 C                 -1.11207400   -1.90010500    0.79689600
 C                 -4.46993800   -0.84121200    0.09330600
 H                 -3.78444300   -4.84155900    0.78073800
 H                 -1.37354400   -5.27811600    1.22660800
 H                  0.27896200   -3.47401600    1.23415700
 N                 -0.57169800    0.47946900    0.57181400
 C                 -6.87460400   -0.53113800   -0.40655800
 C                 -5.55005100    0.05012800   -0.23538600
 C                 -5.39636000    1.41979800   -0.39948300
 C                 -6.48692500    2.23163100   -0.73798400
 C                 -7.77486200    1.64328300   -0.89972800
 C                 -7.99786100    0.29861100   -0.74352800
 H                 -4.42354700    1.87569000   -0.26865500
 H                 -8.60812100    2.28741800   -1.15646000
 C                 -6.30702600    3.69662900   -0.96464900
 H                 -5.41853200    4.07210800   -0.45647300
 H                 -6.18499900    3.88616800   -2.03859000
 H                 -7.18468500    4.25744000   -0.63768600
 C                 -9.34417200   -0.33326200   -0.90491400
 H                 -9.32915600   -1.08128300   -1.70266100
 H                 -9.63888500   -0.85280700    0.01096900
 H                -10.09220300    0.42269400   -1.13994800
 O                 -7.07992400   -1.80715100   -0.25966900
 H                 -6.16498100   -2.26257800   -0.01787900
 N                 -3.15278500   -0.51632400    0.26616900
 N                 -4.66732700   -2.15937100    0.25339000
 H                 -2.66529900    0.37707400    0.21540500
 N                  1.15720800   -0.87579800    1.02613600
 H                  1.73125600   -1.69266400    1.20431600
 C                  4.06630300   -0.05848200    1.42033300
 N                  3.92242400   -1.38223300    1.46163500
 N                  5.35477500    0.32940600    1.69612700
 C                  7.43509800   -0.98299900    2.26936600
 C                  7.85313500   -2.30818200    2.43722400
 C                  6.96808900   -3.38971100    2.28075000
 C                  5.62729600   -3.19696500    1.95149900
 C                  5.17797300   -1.88136700    1.77654600
 C                  6.08489000   -0.80795200    1.93080700
 H                  8.88809700   -2.49843100    2.70062300
 H                  7.34082300   -4.39714200    2.42285800
 H                  4.94878900   -4.03425300    1.84031600
 H                  5.76227900    1.25439100    1.77681100
 C                  8.33709000    0.15123000    2.46061300
 H                  9.37124500   -0.10619900    2.70412700
 N                  7.92244000    1.35921400    2.35987600
 C                  8.84357200    2.47680400    2.56895700
 C                 10.29514300    2.13329300    2.92133400
 C                  8.78653600    3.39178100    1.33416600
 H                  8.42147100    3.04863900    3.40834600
 C                 11.12303300    3.40519700    3.13792000
 H                 10.73408500    1.54782000    2.10413800
 H                 10.32812300    1.50560000    3.81728900
 C                  9.61764500    4.66182800    1.54240300
 H                  9.17171200    2.83199500    0.47285700
 H                  7.74255600    3.63820200    1.11797100
 C                 11.06622500    4.32934000    1.91686000
 H                 12.15941800    3.13769100    3.36625500
 H                 10.73637400    3.94205200    4.01330500
 H                  9.58972100    5.27975000    0.63974300
 H                  9.16708900    5.25597800    2.34719300
 H                 11.62608200    5.24908700    2.11346200
 H                 11.55341700    3.83280700    1.06835200

radii=bondi
dis
rep
cav



