%chk=cyclohexyl-monoBIP-TS2.chk
%mem=10GB
%nprocshared=16
#p opt=(calcfc,tight,ts,noeigentest) freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

cyclohexylimine monoBIP oxidized, TS opt to both protons transferred

1 2
 C                 -0.56914400    2.97565400    0.35301000
 C                 -1.96117400    3.10401200    0.27712700
 C                 -2.54275200    4.37257700    0.19122800
 C                 -1.67015900    5.46371100    0.18645100
 C                 -0.26876500    5.32385800    0.25729800
 C                  0.31740500    4.05965400    0.34088200
 C                 -1.37090400    0.94509300    0.37630400
 H                 -3.61431800    4.50865100    0.12738000
 H                 -2.08765000    6.46097100    0.12232800
 H                  0.35840000    6.20711200    0.24021400
 N                 -2.42562800    1.80208200    0.29674700
 C                 -2.82004000   -1.11126400    0.29828100
 C                 -1.49327700   -0.49288400    0.40706500
 C                 -0.36585200   -1.28607700    0.53730100
 C                 -0.45894700   -2.68984600    0.56398600
 C                 -1.73203400   -3.30120900    0.44695400
 C                 -2.88921100   -2.57194400    0.31446600
 H                  0.61020800   -0.82322400    0.62404900
 H                 -1.79090600   -4.38397600    0.47095900
 C                  0.76672000   -3.53671500    0.70029400
 C                 -4.23702100   -3.20653500    0.19623800
 O                 -3.87013000   -0.42832500    0.20095300
 N                 -0.22418100    1.65191100    0.41490700
 C                  1.73862800    3.75062500    0.38515200
 H                  2.46108400    4.56566000    0.36887000
 N                  2.10297100    2.51277500    0.43436700
 C                  3.52267200    2.13259900    0.44053600
 C                  4.18449600    2.46606900   -0.90678400
 C                  4.29448500    2.77059500    1.60582400
 H                  3.53724500    1.04683400    0.57203500
 C                  5.65160500    2.02566400   -0.91483000
 H                  4.11956200    3.54976300   -1.06233200
 H                  3.62718800    1.98915200   -1.71775400
 C                  5.75991600    2.32446300    1.59275900
 H                  4.24618100    3.86140400    1.51003900
 H                  3.81309600    2.50502400    2.55082700
 C                  6.42964800    2.64394500    0.25168600
 H                  6.11308400    2.29437500   -1.86937400
 H                  5.69828700    0.93253300   -0.83746900
 H                  6.29882900    2.80285900    2.41580200
 H                  5.80775600    1.24298000    1.77037100
 H                  7.46194300    2.28168600    0.24963700
 H                  6.47289400    3.73202100    0.12027500
 H                 -3.37477000    1.43620600    0.24490800
 H                  0.99346700    1.67509100    0.43861100
 C                  1.99587323   -2.69948029    1.10006940
 H                  2.26963108   -2.05790798    0.28868687
 H                  1.75975264   -2.10702105    1.95921963
 H                  2.81225770   -3.35179711    1.33006762
 C                  0.55082739   -4.64614094    1.74631717
 H                 -0.10106493   -5.39337292    1.34434605
 H                  1.49210623   -5.09021062    1.99470449
 H                  0.11143761   -4.22567383    2.62668301
 C                  1.07817640   -4.24004221   -0.63384458
 H                  2.03693025   -4.71124031   -0.57331841
 H                  0.32954044   -4.97869466   -0.83091501
 H                  1.08445995   -3.51886683   -1.42426635
 C                 -5.17884520   -2.68446656    1.29715204
 H                 -5.76579460   -1.87735538    0.91119615
 H                 -5.82555941   -3.47445412    1.61743377
 H                 -4.59956564   -2.33885446    2.12774668
 C                 -4.13048012   -4.73903977    0.30430364
 H                 -3.50610890   -5.10947112   -0.48172511
 H                 -3.70573411   -5.00118210    1.25075608
 H                 -5.10557227   -5.17125791    0.21896447
 C                 -4.89384743   -2.84342680   -1.14850416
 H                 -5.54607699   -3.63541376   -1.45223638
 H                 -5.45685362   -1.94024794   -1.03808685
 H                 -4.13497788   -2.70232934   -1.88952272

radii=bondi
dis
rep
cav













