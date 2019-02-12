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
ref = ('di','E3PT',0.30)
#ref = ('tri','E4PT',0.24)

print('Using reference midpoint potential from: ','cyclohexyl-'+ref[0]+'-BIP ('+ref[1]+')','('+str(ref[2])+' V)') 

A = ('amide-mono','E0PT')

dGc = (gibbs(ref[0],'neutral') + gibbs(A[0],A[1]) - gibbs(A[0],'neutral') - gibbs(ref[0],ref[1]))

Eref = ref[2]

print('E('+A[0]+'-BIP-'+A[1]+'):  \t','{:.3f} V'.format(dGc*27.2114 + Eref))



