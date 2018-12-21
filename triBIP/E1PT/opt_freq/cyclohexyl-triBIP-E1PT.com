%chk=cyclohexyl-triBIP-E1PT.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine triBIP oxidized, first proton transferred, optimization

1 2
 C                  1.65206700    0.37067200    0.92818800
 C                  0.54483800    1.20262000    0.63952500
 C                  0.73158200    2.58322300    0.49612600
 C                  2.02366000    3.08052200    0.64755000
 C                  3.11439500    2.24275700    0.93870600
 C                  2.95778800    0.86011000    1.08874100
 C                 -0.20122000   -0.82249700    0.78172100
 H                 -0.10237000    3.23727100    0.27388700
 H                  2.19959100    4.14398400    0.53750700
 H                  4.10130200    2.67912000    1.03854700
 C                 -2.47161300   -1.75935300    0.55561000
 C                 -3.39925900   -2.81523800    0.56754900
 C                 -3.00761400   -4.12988800    0.82220400
 C                 -1.65195100   -4.33408100    1.05947000
 C                 -0.71989600   -3.28031800    1.04717000
 C                 -1.10404900   -1.95988800    0.79805700
 C                 -4.48179900   -0.91162700    0.10162100
 H                 -3.72218700   -4.94212000    0.83406600
 H                 -1.29693400   -5.33650400    1.26242900
 H                  0.32198900   -3.50333800    1.23869600
 N                 -0.60096700    0.42755800    0.55527700
 C                 -6.91230600   -0.55222400   -0.38525300
 C                 -5.55663200   -0.00937400   -0.22882100
 C                 -5.32587600    1.34325000   -0.39873200
 C                 -6.37311700    2.22626000   -0.73241900
 C                 -7.68445200    1.71739500   -0.88785600
 C                 -7.98031300    0.38314600   -0.73282200
 H                 -4.32497800    1.74331900   -0.28151300
 H                 -8.48009400    2.41093600   -1.13844400
 C                 -6.10081600    3.68079600   -0.94825400
 C                 -9.35795000   -0.16955000   -0.89817600
 O                 -7.16461200   -1.77246700   -0.22959200
 N                 -3.17604100   -0.60868800    0.26710300
 N                 -4.62360200   -2.23896400    0.28459500
 H                 -2.70617800    0.29457100    0.20869600
 N                  1.14743800   -0.90336200    1.00887600
 H                  1.73557100   -1.71016800    1.18997100
 C                  4.04163500   -0.05711400    1.39846500
 N                  3.90640300   -1.38193300    1.43692900
 N                  5.32583900    0.33862100    1.68272800
 C                  7.41115700   -0.96370000    2.26584500
 C                  7.83557200   -2.28688700    2.43308700
 C                  6.95834200   -3.37346900    2.26778300
 C                  5.61860900   -3.18724700    1.93133100
 C                  5.16298600   -1.87364200    1.75773600
 C                  6.06193100   -0.79477300    1.92003700
 H                  8.87015300   -2.47176100    2.70180700
 H                  7.33615700   -4.37909700    2.40921600
 H                  4.94519100   -4.02764800    1.81340800
 H                  5.72563200    1.26655000    1.76730100
 C                  8.30790800    0.17346800    2.46259600
 H                  9.34191600   -0.08116700    2.70938200
 N                  7.89036700    1.38029000    2.36107400
 C                  8.80983600    2.49920900    2.57131100
 C                 10.26135200    2.15692000    2.92535800
 C                  8.75333700    3.41281900    1.33546000
 H                  8.38620700    3.07119000    3.40979200
 C                 11.08902800    3.42913800    3.14039100
 H                 10.70100200    1.57032100    2.10935100
 H                 10.29400200    1.53053700    3.82225700
 C                  9.58439500    4.68311900    1.54228600
 H                  9.13895400    2.85194800    0.47506000
 H                  7.70939900    3.65881800    1.11844900
 C                 11.03283400    4.35102600    1.91764900
 H                 12.12526700    3.16189600    3.36972300
 H                 10.70190700    3.96745200    4.01465500
 H                  9.55680000    5.29969800    0.63868600
 H                  9.13363300    5.27842700    2.34609500
 H                 11.59294900    5.27095800    2.11258500
 H                 11.52002900    3.85269700    1.07017700
 H                 -5.54992500   -2.65783400    0.19363800
 C                 -9.39556935   -1.19383936   -2.04753064
 H                 -9.30157998   -2.18195162   -1.64790561
 H                -10.32448313   -1.10791880   -2.57158870
 H                 -8.58678302   -1.00332892   -2.72167487
 C                -10.36908043    0.95402344   -1.19279095
 H                -11.35088582    0.53713672   -1.27742927
 H                -10.35199391    1.66733593   -0.39542445
 H                -10.10689993    1.43826201   -2.11021909
 C                 -9.79416957   -0.92331717    0.37192146
 H                 -9.50497815   -1.95055018    0.29407122
 H                 -9.32361405   -0.48410496    1.22665542
 H                -10.85700403   -0.85901747    0.47750766
 C                 -5.92176892    3.95971190   -2.45216492
 H                 -6.88117824    4.09187411   -2.90710109
 H                 -5.33944557    4.84773500   -2.58336892
 H                 -5.42028014    3.13331867   -2.91094945
 C                 -7.28324897    4.54011344   -0.46345939
 H                 -6.97921878    5.56430685   -0.40435348
 H                 -8.09670371    4.44873424   -1.15254663
 H                 -7.59538464    4.20435718    0.50335965
 C                 -4.81245366    4.12402986   -0.23043353
 H                 -4.82685657    3.77000910    0.77920099
 H                 -3.96310900    3.71771061   -0.73878214
 H                 -4.75223495    5.19233096   -0.23297425

radii=bondi
dis
rep
cav







