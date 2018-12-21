import cclib
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os
from tqdm import tqdm

class Mode(object):
    def __init__(self,number):
       self.number = number
       self.internal = [] 

class InternalMode(object):
    def __init__(self,name,definition,value,relweight):
        self.name = name
        self.definition = definition
        self.value = value
        self.relweight = relweight

class IRSpec(object):
    def __init__(self,logfile):
        self.data  = cclib.io.ccread(logfile)
        self.scale = 0.9631 # fit to triBIP protonated imine 
        self.freqs = self.data.vibfreqs*self.scale
        self.irs   = self.data.vibirs
        self.modes = dict()
        self.FWHM  = 12
        self.xlim  = (0,4000,0.5)
        self.curvefit()
        self.getNormalModes(logfile)

    def freq(self,number):
        return self.modes[number][0]

    def ir(self,number):
        return self.modes[number][1]

    def components(self,number):
        self.modes[number][2].internal.sort(key=lambda x: x.relweight, reverse=True)
        return self.modes[number][2].internal

    def report(self,window=[0,4000],irthresh=0.0,num_components=5,outfile='output.txt'):
        with open(outfile,'w') as f:
            for index,mode in enumerate(self.modes):
                idx = index+1
                # dont analyze weak peaks
                if self.ir(idx) < irthresh:
                    continue
                # analyze within range 
                elif (self.freq(idx) < min(window)) or (self.freq(idx) > max(window)):
                    continue
                else:
                    self.modes[idx][2].internal.sort(key=lambda x: x.relweight, reverse=True)
                    f.write('Mode %3d: %10.2f %6.1f\n' % (idx,self.freq(idx),self.ir(idx)) )
                    f.write('%10s %16s %10s %20s\n' % ('Name','Definition','Value','Relative Weight (%)'))
                    # check that we don't have fewer than desired components
                    if len(self.modes[idx][2].internal) < num_components: 
                        length = len(self.modes[idx][2].internal)
                    else:
                        length = num_components
                    for i in range(length):
                        component = self.modes[idx][2].internal[i]
                        f.write('%10s %16s %10.2f %20.2f\n' % (component.name,component.definition,component.value,component.relweight))


    def topcomponents(self,number):
        self.modes[number][2].internal.sort(key=lambda x: x.relweight, reverse=True)
        for i in range(4):
            component = self.modes[number][2].internal[i]
            print(component.name,component.definition,component.value,component.relweight)

    def lorentzian(self,x0,h,x,gm):
         return (h*np.power(gm/2,2))/(np.power(x-x0,2) + np.power(gm/2,2))

    def curvefit(self):
        bands = zip(self.freqs,self.irs)
        self.x = np.arange(self.xlim[0],self.xlim[1],self.xlim[2])
        self.curve = np.zeros_like(self.x)
        for band in bands:
            self.curve += self.lorentzian(band[0],band[1],self.x,self.FWHM)


    def getNormalModes(self,logfile):
        with open(logfile,'r') as f:
            copy = False
            number = None
            for line in f.readlines():
                try:
                    a,b,c = line.split()[1:4]
                except (IndexError, ValueError):
                    a,b,c = None,None,None 
                if ('Normal','Mode') == (a,b):
                        number = int(c)
                        copy = True
                        self.modes[number] = [self.freqs[number-1],self.irs[number-1],Mode(number)]
                elif ['Center', 'Atomic', 'Forces', '(Hartrees/Bohr)'] == line.split():
                    copy = False
                elif copy and (line.split()[0] == '!') \
                        and (line.split()[1][0] in ['A','D','R']) \
                        and (line.split()[-1] == '!'):
                    name = line.split()[1]
                    definition = line.split()[2]
                    value = float(line.split()[3])
                    relweight = float(line.split()[4])
                    self.modes[number][2].internal.append(InternalMode(name,definition,value,relweight)) 

            
if __name__ == '__main__':

    prefixes = ['mono','di','tri']
    states   = ['neutral','E0PT','E1PT','E2PT','E3PT','E4PT']
    # tqdm is just the progress bar
    for prefix, state in tqdm(itertools.product(prefixes,states),total=len(prefixes)*len(states)): 
        logfile = '../../'+prefix+'BIP/'+state+'/opt_freq/'+'cyclohexyl-'+prefix+'BIP-'+state+'.log'
        if os.path.isfile(logfile):
            mol = IRSpec(logfile)
            mol.report(window=[1000,4000],irthresh=20.0,outfile=prefix+'BIP/cyclohexyl-'+prefix+'BIP-'+state+'-modes.txt')
            #plt.plot(mol.x,mol.curve)    
            #plt.xlim([3700,1000])
            #plt.show()


