%chk=unsub-monoBIP-E0PT.chk
%mem=10GB
%nprocshared=16
#p opt=tight freq=internalmodes b3lyp/6-31g(d,p)
scrf=(cpcm,solvent=dichloromethane,read) nosymm empiricaldispersion=gd3bj
int=grid=ultrafine

Simple unsub-monoBIP E0PT optimization

1 2
 C                 -0.48305300    2.71558400    0.41614500
 C                 -1.88105200    2.89835900    0.32760100
 C                 -2.41289400    4.19012200    0.24342100
 C                 -1.52045600    5.26198800    0.25627700
 C                 -0.13359500    5.06453100    0.35031400
 C                  0.42705100    3.78268800    0.42980400
 C                 -1.51854300    0.75201500    0.42465100
 H                 -3.48207900    4.34633900    0.16902100
 H                 -1.90526600    6.27315300    0.19413200
 H                  0.52953900    5.92276100    0.36061700
 N                 -2.49642700    1.65628200    0.33630200
 C                 -3.06366800   -1.18903700    0.35231400
 C                 -1.73921900   -0.68500300    0.44726500
 C                 -0.65872400   -1.57362300    0.54049700
 C                 -0.84875700   -2.95028900    0.54135500
 C                 -2.16787400   -3.41764100    0.44644700
 C                 -3.28612600   -2.58500100    0.35332600
 H                  0.34704400   -1.17537400    0.61268000
 H                 -2.33231200   -4.48382500    0.44549100
 C                  0.36142000   -3.88779800    0.61444200
 C                 -4.70520000   -3.16741300    0.25338500
 O                 -4.11725500   -0.34687600    0.25749700
 H                 -3.76780900    0.59151100    0.26864100
 N                 -0.28701400    1.35660000    0.47643900
 H                  0.63826000    0.94763200    0.52815200
 C                  1.24577700   -3.65452100   -0.62728500
 H                  1.57556200   -2.61430600   -0.68350800
 H                  2.13381600   -4.29237600   -0.58864300
 H                  0.69738400   -3.88708200   -1.54417900
 C                  1.18278600   -3.57554900    1.88073100
 H                  0.58865300   -3.74388700    2.78286100
 H                  2.06432100   -4.22163400    1.92638700
 H                  1.52489900   -2.53762700    1.88737400
 C                 -0.04946000   -5.36608800    0.65210000
 H                 -0.66468200   -5.59166500    1.52739600
 H                 -0.60816100   -5.65207000   -0.24333100
 H                  0.84524000   -5.99163900    0.70128300
 C                 -5.35900200   -2.73216100   -1.07525100
 H                 -6.37293400   -3.13822900   -1.13995800
 H                 -5.41369800   -1.64726500   -1.15535500
 H                 -4.78635400   -3.11381600   -1.92527500
 C                 -5.55747400   -2.68254900    1.44568700
 H                 -5.65813200   -1.59820600    1.44662500
 H                 -6.55654600   -3.12576700    1.39415100
 H                 -5.09795700   -2.98687900    2.39027700
 C                 -4.69423100   -4.70465900    0.28381600
 H                 -4.12605200   -5.12340000   -0.55105600
 H                 -4.26883500   -5.08975800    1.21438200
 H                 -5.72145900   -5.06913900    0.20754000
 H                  1.48364318    3.62766802    0.49675010

radii=bondi
dis
rep
cav









