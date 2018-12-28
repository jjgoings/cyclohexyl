import cclib
#import seaborn as sns
import numpy as np
#import matplotlib.pyplot as plt
import os

#sns.set_color_codes()

class Structure(object):
    def __init__(self,logfile):
        self.data     = cclib.io.ccread(logfile)
        self.charge   = self.data.charge
        self.mult     = self.data.mult
        self.vibfreqs = self.data.vibfreqs
        self.natom   = self.data.natom
        self.pressure = self.data.pressure
        self.temperature = self.data.temperature
        self.atommasses = self.data.atommasses
        self.nelectrons = self.data.nelectrons
        if 'opt_freq' in logfile:
            self.optdone  = self.data.optdone
    
class Molecule(object):
    def __init__(self,prefix):
        self.prefix  = str(prefix)
        path = '../../'+self.prefix+'BIP/'
        self.states = []

        for state in ['neutral','E0PT','E1PT','E2PT','E3PT','E4PT']:
            logfile = path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'
            if os.path.isfile(logfile):
                setattr(self,state.lower(),Structure(path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'))
                setattr(getattr(self,state.lower()),'thermo',Structure(path+state+'/thermo/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'))
                setattr(getattr(self,state.lower()),'name',self.prefix+'BIP-cyclohexylimine '+'('+state+')')
                self.states.append(state)
               
    def check_structure(self):
        print('Checking ',self.prefix+'BIP-cyclohexylimine structures...')
        for state,data in [(x,getattr(self,x.lower())) for x in self.states]: 
            print(' > '+state)
            # Correct charge/multiplicity?
            if state == 'neutral':
                assert data.charge == data.thermo.charge == 0
                assert data.mult == data.thermo.mult  == 1
            else:
                assert data.charge == data.thermo.charge == 1
                assert data.mult == data.thermo.mult  == 2
            print("   [+] Charge and multiplicity...OK")
            

            assert data.optdone # Converged?
            assert all(v > 0 for v in [min(data.vibfreqs), min(data.thermo.vibfreqs)]) # Local minimum? 
            print("   [+] Converged to local minima...OK")

            # thermochem
            assert (data.temperature == data.thermo.temperature == 298.15) # K
            assert (data.pressure == data.thermo.pressure == 1.0) # atm
            print("   [+] Thermochemistry at 298.15K and 1atm...OK")

            # are frequencies consistent?
            assert np.allclose(data.vibfreqs,data.thermo.vibfreqs) # Local minimum? 
            print("   [+] Opt and Thermo have same unscaled frequencies?...OK")
 

            if self.prefix == 'mono': 
                assert data.natom == data.thermo.natom == 69 
                # note: masses use 2012 CODATA recommendations
                assert np.allclose(sum(data.atommasses[:data.natom]),431.2936616)
                assert (data.nelectrons + data.charge) == (data.thermo.nelectrons + data.thermo.charge) == 234 
            elif self.prefix == 'di': 
                assert data.natom == data.thermo.natom == 82 
                assert np.allclose(sum(data.atommasses[:data.natom]),547.3311096)
                assert (data.nelectrons + data.charge) == (data.thermo.nelectrons + data.thermo.charge) == 294 
            elif self.prefix == 'tri': 
                assert data.natom == data.thermo.natom == 95 
                assert np.allclose(sum(data.atommasses[:data.natom]),663.3685576)
                assert (data.nelectrons + data.charge) == (data.thermo.nelectrons + data.thermo.charge) == 354 
            else: 
               print("ERROR: No prefix match") 
            # consistent masses?
            assert np.allclose(sum(data.atommasses[:data.natom]),sum(data.atommasses[:data.natom]))
            print("   [+] Correct number of atoms, molecular weights, and electrons...OK")

if __name__ == '__main__':
    mono = Molecule('mono')
    mono.check_structure()
    di = Molecule('di')
    di.check_structure() 
    tri = Molecule('tri')
    tri.check_structure()
    print("Structures look good!")
    


