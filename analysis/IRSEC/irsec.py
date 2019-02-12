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
        self.FWHM  = 12
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
                if state == 'neutral':
                    setattr(getattr(self,state.lower()),'name','['+self.prefix+'-BIP-cyclohexylimine]')
                else:
                    setattr(getattr(self,state.lower()),'name','['+self.prefix+'-BIP-cyclohexylimine]$^{+\\bullet}$')
  
        self.freq  = self.neutral.x # assumes neutral always exists! 

    def plotSpectra(self,curves=[],interpolate=False,num=5,xlim=None,save=None,show=True,colors=None):
        curves = [curves] if isinstance(curves, IRSpec) else curves
        if not xlim:
            xlim = [max(self.freq),min(self.freq)]
        fig,ax = plt.subplots(figsize=(8,6))
        ax.tick_params(width=2)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(2)

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
                    if not colors:
                        color = 'r'
                    else:
                        color = colors
                    zorder = 3
                else:
                    label = None
                    color = 'gray'
                    zorder = 1
                plt.plot(self.freq,curve,label=label,color=color,zorder=zorder,lw=2)
           #plt.ylim([-50,max(max(start.curve),max(finish.curve))])
            plt.ylim([-50,1200])
            #w = 3
            #plt.axvspan(1659+w, 1638-w, facecolor='r', alpha=0.2,zorder=-1,label='imine C=N str')
            #plt.axvspan(1616+w, 1563-w, facecolor='orange', alpha=0.2,zorder=-1,label='aromatic sym C=C str')
            #plt.axvspan(1532+w, 1502-w, facecolor='y', alpha=0.2,zorder=-1,label='interaromatic C-C str')
            #plt.axvspan(1493+w, 1477-w, facecolor='g', alpha=0.2,zorder=-1,label='benzimidazole$_1$ NH bend + phenol C=O str')
            #plt.axvspan(1471+w, 1376-w, facecolor='b', alpha=0.2,zorder=-1,label='t-butyl C-H bend')
            #plt.axvspan(1428+w, 1366-w, facecolor='purple', alpha=0.2,zorder=-1,label='aromatic asym C=C str')
            #plt.axvspan(1355+w, 1334-w, facecolor='pink', alpha=0.2,zorder=-1,label='benzimidazole C=N str')

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
                plt.ylim([-50,1300])
            
            
        plt.legend(fontsize=14)
        plt.xlim(xlim)
        plt.xticks(np.arange(min(xlim),max(xlim),50),fontsize=14)
        ax.set_xticks(np.arange(min(xlim),max(xlim),25),minor=True)
       
        plt.yticks([0],fontsize=14)
        #plt.grid(color='gray',alpha=0.4)
        plt.grid(color='gray',alpha=0.4,which='minor',ls='--')
        plt.xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)',fontsize=24)
        plt.ylabel('Intensity (arb. units)',fontsize=24)
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
    mono.plotSpectra(curves=[mono.neutral,mono.e2pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='monoIRSEC.pdf',colors='b')
    di = Molecule('di')
    di.plotSpectra(curves=[di.neutral,di.e3pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='diIRSEC.pdf',colors='g')
    tri = Molecule('tri')
    tri.plotSpectra(curves=[tri.neutral,tri.e4pt],interpolate=True,num=5,xlim=[1700,1300],show=False,save='triIRSEC.pdf',colors='r')
    mono = Molecule('mono')
    mono.plotSpectra(curves=[mono.neutral,mono.e2pt],interpolate=True,num=5,xlim=[3500,3200],show=False,save='monoIRSEC_3500.pdf',colors='b')
    di = Molecule('di')
    di.plotSpectra(curves=[di.neutral,di.e3pt],interpolate=True,num=5,xlim=[3500,3200],show=False,save='diIRSEC_3500.pdf',colors='g')
    tri = Molecule('tri')
    tri.plotSpectra(curves=[tri.neutral,tri.e4pt],interpolate=True,num=5,xlim=[3500,3150],show=False,save='triIRSEC_3500.pdf',colors='r')


