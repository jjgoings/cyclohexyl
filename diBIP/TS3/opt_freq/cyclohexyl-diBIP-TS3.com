%chk=cyclohexyl-diBIP-TS3.chk
%mem=10GB
%nprocshared=16
#p opt=(calcall,tight,ts,noeigentest) freq=(readfc,internalmodes) b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm guess=read
empiricaldispersion=gd3bj int=grid=ultrafine

cyclohexylimine diBIP oxidized, last transition state optimization

1 2
 C                  1.66887300    1.72684900    0.48892900
 C                  0.43215400    2.36813600    0.36153000
 C                  0.36944900    3.76309500    0.35057400
 C                  1.58160700    4.45160000    0.47272800
 C                  2.82156800    3.79584800    0.59564700
 C                  2.88968000    2.39961200    0.60530100
 C                  0.17933000    0.14169400    0.31840500
 H                 -0.56847500    4.29327500    0.24940300
 H                  1.56443100    5.53446800    0.46654300
 H                  3.73138700    4.37886000    0.67670600
 C                 -1.88743600   -1.20877800    0.09053600
 C                 -2.52181000   -2.47440400   -0.00900100
 C                 -1.80965500   -3.67224500    0.02408900
 C                 -0.42492900   -3.57948800    0.15985100
 C                  0.22908300   -2.34162700    0.25901100
 C                 -0.48213200   -1.13756700    0.22443800
 C                 -3.99843100   -0.82013100   -0.09754000
 H                 -2.30741300   -4.63059400   -0.05277400
 H                  0.16506700   -4.48691200    0.18964400
 H                  1.30746400   -2.31345000    0.36041300
 N                 -0.48288000    1.33071900    0.25777500
 C                 -6.51084100   -0.95093700   -0.36233800
 C                 -5.28487200   -0.16106700   -0.19629000
 C                 -5.36087100    1.22338000   -0.13957100
 C                 -6.58865400    1.89416400   -0.23758300
 C                 -7.78306700    1.13994400   -0.38332400
 C                 -7.78059500   -0.23084300   -0.44435100
 H                 -4.44779800    1.79440000   -0.02193800
 H                 -8.72427700    1.67558300   -0.45070800
 C                 -6.66257200    3.38882400   -0.17261300
 H                 -7.02189100    3.70729400    0.81370900
 H                 -5.68757600    3.84550300   -0.34205100
 H                 -7.37878100    3.77166000   -0.90414800
 C                 -9.03005100   -1.03860200   -0.60192500
 H                 -9.04526000   -1.54569100   -1.57162100
 H                 -9.08017400   -1.82447000    0.15632300
 H                 -9.91483700   -0.40657400   -0.52510500
 O                 -6.49468600   -2.20722200   -0.43762800
 N                 -2.81757100   -0.19989500    0.03097400
 N                 -3.85783500   -2.18096200   -0.12647600
 N                  1.49592500    0.36535400    0.46409800
 C                  4.09158800    1.58051900    0.68780800
 H                  5.06269600    2.06797600    0.74169900
 N                  3.96545600    0.29648700    0.67588900
 C                  5.09305000   -0.63777600    0.69028100
 C                  4.86608500   -1.66219600   -0.43394200
 C                  6.47163200    0.01353800    0.57746900
 H                  5.03077900   -1.17176700    1.64744100
 C                  5.97769000   -2.71480200   -0.45265800
 H                  4.84442800   -1.12415800   -1.38868000
 H                  3.88682400   -2.13294500   -0.30013200
 C                  7.57958800   -1.04632000    0.55273300
 H                  6.50666300    0.60127900   -0.34737100
 H                  6.63486600    0.70738200    1.40735400
 C                  7.35888400   -2.06414400   -0.57093300
 H                  5.80972400   -3.41557900   -1.27530400
 H                  5.93499200   -3.29820900    0.47500200
 H                  8.54967100   -0.55416400    0.43704800
 H                  7.60089300   -1.56889700    1.51681800
 H                  8.14037900   -2.83026100   -0.54544300
 H                  7.43910400   -1.55643100   -1.54027200
 H                 -4.66603200   -2.78935000   -0.22753400
 H                 -1.49835200    1.36186000    0.15903500
 H                  2.63511400   -0.06556500    0.56473600


radii=bondi
dis
rep
cav









