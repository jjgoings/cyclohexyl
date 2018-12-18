%chk=cyclohexyl-triBIP-E4PT.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine triBIP oxidized, all protons transferred optimization

1 2
 C                  1.63303900    0.41917900    0.99157600
 C                  0.56491000    1.30227900    0.69444700
 C                  0.74832300    2.67709600    0.56558200
 C                  2.04585800    3.16112900    0.74306000
 C                  3.11781900    2.30862100    1.04121700
 C                  2.93828500    0.92585700    1.17090400
 C                 -0.13574000   -0.80167200    0.81952800
 H                 -0.07606600    3.34046200    0.33719900
 H                  2.22992800    4.22456700    0.64787300
 H                  4.11280700    2.71530000    1.17146100
 C                 -2.44704000   -1.70701400    0.53638700
 C                 -3.33861100   -2.81289500    0.53761800
 C                 -2.91602900   -4.11926400    0.78325400
 C                 -1.55769400   -4.29961100    1.03326500
 C                 -0.65454100   -3.22304800    1.03962700
 C                 -1.07230900   -1.91159200    0.79870800
 C                 -4.40543500   -0.91918800    0.08812800
 H                 -3.61005000   -4.95038900    0.78116700
 H                 -1.18213900   -5.29671900    1.22933000
 H                  0.39443800   -3.40246200    1.23837400
 N                 -0.54168900    0.48694300    0.59535600
 C                 -6.86085700   -0.58200000   -0.41215900
 C                 -5.50788800   -0.03830800   -0.24080700
 C                 -5.29072100    1.32311800   -0.39654900
 C                 -6.33502500    2.19968200   -0.72706300
 C                 -7.64684800    1.68530000   -0.89592100
 C                 -7.93262300    0.35024000   -0.75428100
 H                 -4.28807600    1.71127300   -0.26172100
 H                 -8.44566700    2.37704300   -1.14374000
 C                 -6.07392300    3.66069000   -0.92898400
 H                 -5.15594600    3.97220800   -0.42886700
 H                 -5.96383400    3.87877800   -1.99791000
 H                 -6.90797200    4.26321500   -0.56253900
 C                 -9.30897800   -0.20937900   -0.92757200
 H                 -9.33271700   -0.93476300   -1.74667300
 H                 -9.62307000   -0.74742900   -0.02868200
 H                -10.02698400    0.58356800   -1.13646200
 O                 -7.11727600   -1.80718700   -0.27375900
 N                 -3.13048800   -0.54528200    0.25567300
 N                 -4.56573100   -2.26875300    0.25294100
 N                  1.17343400   -0.88094400    1.06044300
 C                  4.04222700    0.03648200    1.47486500
 N                  3.86989900   -1.31839900    1.54049400
 N                  5.31364100    0.39243500    1.70684300
 C                  7.34243300   -1.00829900    2.22665700
 C                  7.76617600   -2.33990300    2.39910300
 C                  6.87119600   -3.40690800    2.29077700
 C                  5.51623000   -3.20033200    2.00340100
 C                  5.09127000   -1.88779400    1.82532500
 C                  5.97845300   -0.78786300    1.92472500
 H                  8.80785200   -2.53159100    2.62907200
 H                  7.23245800   -4.41736900    2.43349000
 H                  4.82724700   -4.03228800    1.92675300
 C                  8.26257900    0.07310500    2.39903400
 H                  9.29603300   -0.14935900    2.63661200
 N                  7.91396300    1.32092000    2.30809300
 C                  8.77405300    2.49725400    2.51550900
 C                 10.20784100    2.13165400    2.88900000
 C                  8.72140500    3.39085300    1.26896800
 H                  8.31966800    3.04410600    3.35009800
 C                 11.03652700    3.40332200    3.10851300
 H                 10.65203000    1.54256900    2.07846400
 H                 10.21998800    1.51428600    3.79164900
 C                  9.55082400    4.66010700    1.49316000
 H                  9.12001700    2.82464500    0.41989800
 H                  7.68167200    3.64208400    1.03909700
 C                 10.99398200    4.32482800    1.88495400
 H                 12.06723200    3.12484200    3.34397100
 H                 10.64572300    3.94036500    3.98090400
 H                  9.53208300    5.27403800    0.58891300
 H                  9.08881600    5.25390400    2.29081500
 H                 11.54973800    5.24420000    2.09090200
 H                 11.49225000    3.82940700    1.04287500
 H                  6.92745500    1.49483500    2.07208100
 H                  2.95617800   -1.75331500    1.40587600
 H                 -1.51061400    0.72627500    0.39419700
 H                 -5.47913200   -2.70290000    0.15106500

radii=bondi
dis
rep
cav





