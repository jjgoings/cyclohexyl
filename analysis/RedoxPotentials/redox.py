import cclib
import subprocess 
import os
import numpy as np

def gibbs(prefix,configuration,mode=1):
    if mode == 1:
        path = '../../'+prefix+'BIP/'+configuration+'/opt_freq/'
        filename = path+'cyclohexyl-'+prefix+'BIP-'+configuration+'.log'

        command = ["python", "-m", "goodvibes", 
                    "-q", "grimme",
                    "-v", "0.9649",
                    "-f", "100.0",
                    filename]   
 
        data = subprocess.check_output(command).decode("utf-8")
 
        data = data.split('\n')[-3].split()
        assert((len(data) == 9) and (data[1] == 'cyclohexyl-'+prefix+'BIP-'+configuration))

        freeenergy = float(data[-1])

        return freeenergy
    elif mode == 2:
        path = '../../'+prefix+'BIP/'+configuration+'/thermo/'
        filename = path+'cyclohexyl-'+prefix+'BIP-'+configuration+'.log'
        data = cclib.io.ccread(filename)
        return data.freeenergy

mode = 1
#ref = ('mono','E2PT',0.36)
#ref = ('di','E3PT',0.30)
ref = ('tri','E4PT',0.24)

print('Using reference midpoint potential from: ','cyclohexyl-'+ref[0]+'-BIP','('+str(ref[2])+' V)') 

mono = ('mono','E2PT')
di   = ('di','E3PT')
tri  = ('tri','E4PT')

for A in [mono,di,tri]:
    dGc = (gibbs(ref[0],'neutral',mode) + gibbs(A[0],A[1],mode) - gibbs(A[0],'neutral',mode) - gibbs(ref[0],ref[1],mode))
    Eref = ref[2] 
    print('E(cyclohexyl-'+A[0]+'-BIP-'+A[1]+'):  \t','{:.2f} V'.format(dGc*27.2114 + Eref))


