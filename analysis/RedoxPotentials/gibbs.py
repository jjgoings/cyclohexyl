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
                    "-v", "0.9631",
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

mode = 2

#mono = ('mono','E2PT')
#di   = ('di','E3PT')
#tri  = ('tri','E4PT')

mono0 = ('mono','E0PT')
mono1 = ('mono','E1PT')
mono2 = ('mono','E2PT')
di0   = ('di','E0PT')
di1   = ('di','E1PT')
di2   = ('di','E2PT')
di3   = ('di','E3PT')
tri0  = ('tri','E0PT')
tri1  = ('tri','E1PT')
tri2  = ('tri','E2PT')
tri3  = ('tri','E3PT')
tri4  = ('tri','E4PT')

for A in [tri0,tri1,tri2,tri3,tri4]:
    dGc = gibbs(A[0],A[1],mode) - gibbs('tri','neutral',mode)
    print('E(cyclohexyl-'+A[0]+'-BIP-'+A[1]+'):  \t','{:.2f} kcal/mol'.format((dGc*627.509)))



