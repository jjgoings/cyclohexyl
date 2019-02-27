import cclib
import subprocess 
import os
import numpy as np

def gibbs(prefix,configuration):
    if prefix =='amide-mono':
        path = '../../'+prefix+'BIP/'+configuration+'/thermo/'
        filename = path+prefix+'BIP-'+configuration+'.log'
        data = cclib.io.ccread(filename)
        return data.freeenergy
    else:
        path = '../../'+prefix+'BIP/'+configuration+'/thermo/'
        filename = path+'cyclohexyl-'+prefix+'BIP-'+configuration+'.log'
        data = cclib.io.ccread(filename)
        return data.freeenergy

mode = 2
#ref = ('mono','E2PT',0.36)
ref = ('di','E3PT',0.30 + 0.46) # V vs SCE in DCM, instead of Fc
#ref = ('tri','E4PT',0.24)

print('Using reference midpoint potential from: ','cyclohexyl-'+ref[0]+'-BIP ('+ref[1]+')','('+str(ref[2])+' V)') 

a = ('amide-mono','E0PT')
b = ('amide-mono','E1PT')
c = ('amide-mono','E2PT')
d = ('amide-mono','E2PTz')

for A in [a,b,c,d]:
    dGc = (gibbs(ref[0],'neutral') + gibbs(A[0],A[1]) - gibbs(A[0],'neutral') - gibbs(ref[0],ref[1]))
    
    Eref = ref[2] # V vs SCE in DCM
    
    print('E('+A[0]+'-BIP-'+A[1]+'):  \t','{:.3f} V'.format(dGc*27.2114 + Eref))



