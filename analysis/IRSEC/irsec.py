import cclib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os

sns.set_color_codes()

class IRSpec(object):
    def __init__(self,logfile):
        self.data  = cclib.io.ccread(logfile)
        self.scale = 0.9631 # fit to tri-BIP E4PT protonated imine 
        self.freqs = self.data.vibfreqs*self.scale
        self.irs   = self.data.vibirs
        self.FWHM  = 6
        self.xlim  = (0,4000,0.5)
        self.curvefit()
    
    def lorentzian(self,x0,h,x,gm):
         return (h*np.power(gm/2,2))/(np.power(x-x0,2) + np.power(gm/2,2))
        
    
    def curvefit(self):
        bands = zip(self.freqs,self.irs)
        self.x = np.arange(self.xlim[0],self.xlim[1],self.xlim[2])
        self.curve = np.zeros_like(self.x) 
        for band in bands:
            self.curve += self.lorentzian(band[0],band[1],self.x,self.FWHM) 
     

class Molecule(object):
    def __init__(self,prefix):
        self.prefix  = str(prefix)
        path = '../../'+self.prefix+'BIP/'

        for state in ['neutral','E0PT','E1PT','E2PT','E3PT','E4PT']:
            logfile = path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'
            if os.path.isfile(logfile):
                setattr(self,state.lower(),IRSpec(path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'))
                setattr(getattr(self,state.lower()),'name',self.prefix+'BIP-cyclohexylimine '+'('+state+')')
  
        self.freq  = self.neutral.x # assumes neutral always exists! 

    def plotSpectra(self,curves=[],interpolate=False,num=5,xlim=None,save=None,show=True,colors=None):
        curves = [curves] if isinstance(curves, IRSpec) else curves
        if not xlim:
            xlim = [max(self.freq),min(self.freq)]
        fig,ax = plt.subplots(figsize=(8,6))

        if interpolate:
            start = curves[0]
            finish = curves[-1]
            #ax.set_prop_cycle('color',plt.cm.Reds(np.linspace(0.3,1,num)))
            for step in np.linspace(0,num,num=num):
                curve = ((num - step)/num)*start.curve + (step/num)*finish.curve
                if step == 0:
                    label = start.name
                    color = 'black'
                    zorder = 2
                elif step == np.linspace(0,num,num=num)[-1]:
                    label = finish.name
                    color = 'r'
                    zorder = 3
                else:
                    label = None
                    color = 'gray'
                    zorder = 1
                plt.plot(self.freq,curve,label=label,color=color,zorder=zorder)
            plt.ylim([-50,max(max(start.curve),max(finish.curve))])
            plt.axvspan(1700, 1615, facecolor='r', alpha=0.2,zorder=-1,label='C=N')
            plt.axvspan(1615, 1580, facecolor='orange', alpha=0.2,zorder=-1,label='arom C=C')
            plt.axvspan(1225,  950, facecolor='y', alpha=0.2,zorder=-1,label='arom C-H in-plane bend')
            plt.axvspan(1410, 1310, facecolor='g', alpha=0.2,zorder=-1,label='phenol OH bend')
            plt.axvspan(1210, 1190, facecolor='b', alpha=0.2,zorder=-1,label='phenol C-O stretch')

        # e.g. no interpolation! finish = None
        else:
            if colors: assert len(colors) == len(curves)
            #ax.set_prop_cycle('color',plt.cm.tab10(np.linspace(0.2,0.8,len(curves))))
            #ax.set_prop_cycle('color',plt.cm.tab10(range(10)))
            for idx,spec in enumerate(curves):
                if colors:
                    color = colors[idx]
                else: color = None
                plt.plot(self.freq,spec.curve,label=spec.name,color=color)
                plt.ylim([-50,2000])
            
            
        plt.legend()
        plt.xlim(xlim)
        plt.xticks(np.arange(min(xlim),max(xlim),50))
        ax.set_xticks(np.arange(min(xlim),max(xlim),25),minor=True)
       
        plt.yticks([0])
        #plt.grid(color='gray',alpha=0.4)
        plt.grid(color='gray',alpha=0.4,which='minor',ls='--')
        plt.xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)')
        plt.ylabel('Intensity (arb. units)')
        if save:
            if isinstance(save,str):
                plt.tight_layout()
                plt.savefig(save)
            else:
                raise ValueError("'save' needs to be your filename!")
        if show:
            plt.show()
        plt.close()

    def plotDifference(self,mol1,mol2,xlim=None,save=None,show=True):
        if not xlim:
            xlim = [max(mol1.x),min(mol1.x)]
        fig,ax = plt.subplots()
        label = "$\Delta$ Abs. "+str(mol1.name)+'$-$'+str(mol2.name)
        plt.plot(self.freq,mol1.curve-mol2.curve,label=label)
        plt.legend()
        plt.xlim(xlim)
        plt.xticks(np.arange(min(xlim),max(xlim),500))
        ax.set_xticks(np.arange(min(xlim),max(xlim),100),minor=True)
        #plt.ylim([-50,max(max(start.curve),max(finish.curve))])
        plt.yticks([0])
        plt.grid(color='gray',alpha=0.4)
        plt.grid(color='gray',alpha=0.4,which='minor',ls='--')
        plt.xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)')
        plt.ylabel('$\Delta$ Intensity (arb. units)')
        if save:
            if isinstance(save,str):
                plt.tight_layout()
                plt.savefig(save)
            else:
                raise ValueError("'save' needs to be your filename!")
        if show:
            plt.show()
        plt.close()

    
    def dumpVibs(self,which,thresh=75,output=None):
        """ Return some number peaks from the frequency calculation
            Variables:
            thresh = float. minimum intensity.

            Returns: Energy (cm-1), Intensity 
        """

        idx = np.where(which.irs > thresh)
        peaks = np.hstack((which.freqs[idx].reshape(-1,1),which.irs[idx].reshape(-1,1)))[::-1]
        if output:
            np.savetxt(output,peaks,fmt='%.1f',delimiter=',')

if __name__ == '__main__':
#    tri = Molecule('tri')
#    tri.plotSpectra(curves=[tri.neutral,tri.e4pt],interpolate=True,num=5,xlim=[1700,1300],show=True)

    mono = Molecule('mono')
    mono.plotSpectra(curves=[mono.neutral,mono.e2pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='monoIRSEC.pdf')
    di = Molecule('di')
    di.plotSpectra(curves=[di.neutral,di.e3pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='diIRSEC.pdf')
    tri = Molecule('tri')
    tri.plotSpectra(curves=[tri.neutral,tri.e4pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='triIRSEC.pdf')
    mono = Molecule('mono')
    mono.plotSpectra(curves=[mono.neutral,mono.e2pt],interpolate=True,num=5,xlim=[3700,3000],show=False,save='monoIRSEC_3500.pdf')
    di = Molecule('di')
    di.plotSpectra(curves=[di.neutral,di.e3pt],interpolate=True,num=5,xlim=[3700,3000],show=False,save='diIRSEC_3500.pdf')
    tri = Molecule('tri')
    tri.plotSpectra(curves=[tri.neutral,tri.e4pt],interpolate=True,num=5,xlim=[3700,3000],show=False,save='triIRSEC_3500.pdf')


