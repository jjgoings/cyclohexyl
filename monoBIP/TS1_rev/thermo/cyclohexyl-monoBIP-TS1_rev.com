%chk=cyclohexyl-monoBIP-TS1_rev.chk
%mem=10GB
%nprocshared=16
#p freq(ReadFC,InternalModes) Geom(AllCheck) NoSymm scale=0.9631

