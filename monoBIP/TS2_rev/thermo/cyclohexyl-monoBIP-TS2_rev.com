%chk=cyclohexyl-monoBIP-TS2_rev.chk
%mem=10GB
%nprocshared=4
#p freq(ReadFC,InternalModes) Geom(AllCheck) NoSymm scale=0.9649

