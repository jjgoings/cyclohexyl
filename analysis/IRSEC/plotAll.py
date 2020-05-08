import cclib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
import matplotlib.gridspec as gridspec

matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams.update({'font.size': 18})

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

        if prefix == 'mono': num = str(2)
        if prefix == 'di':   num = str(3)
        if prefix == 'tri':  num = str(4)

        for state in ['neutral','E0PT','E1PT','E2PT','E3PT','E4PT']:
            logfile = path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'
            if os.path.isfile(logfile):
                setattr(self,state.lower(),IRSpec(path+state+'/opt_freq/cyclohexyl-'+self.prefix+'BIP-'+state+'.log'))
                if state == 'neutral':
                    #setattr(getattr(self,state.lower()),'name','['+self.prefix+'-BIP-cyclohexylimine]')
                    setattr(getattr(self,state.lower()),'name','[$\\bf{'+num+'}$]')
                else:
                    #setattr(getattr(self,state.lower()),'name','['+self.prefix+'-BIP-cyclohexylimine]$^{+\\bullet}$')
                    setattr(getattr(self,state.lower()),'name','[$\\bf{'+num+'}$]$^{+\\bullet}$')
  
        self.freq  = self.neutral.x # assumes neutral always exists! 

#    def plotAll(self,curves=[],interpolate=False,num=5,xlim=None,save=None,show=True,colors=None):

interpolate = True
num = 5
show = True
save = 'allIRSEC.png' 

fig, ax = plt.subplots(nrows=3, ncols=2,sharex='col',sharey='row',figsize=(14,8))
#fig.tight_layout()
fig.subplots_adjust(hspace=0,wspace=0.05)
big = fig.add_subplot(111)
#
big.set_zorder(-100)
big.set_facecolor('none')
big.spines['top'].set_color('none')
big.spines['bottom'].set_color('none')
big.spines['left'].set_color('none')
big.spines['right'].set_color('none')
big.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
big.set_xlabel(r'$\tilde{\nu}$ (cm$^{-1}$)',fontsize=32)
big.set_ylabel('Intensity (arb. units)',fontsize=32)




for row,mol in enumerate([Molecule('mono'),Molecule('di'),Molecule('tri')]):
    for col,region in enumerate(['high','low']):
 
        if mol.prefix == 'mono': 
            curves = [mol.neutral,mol.e2pt]
            colors='#2200F5'
        if mol.prefix == 'di':   
            curves = [mol.neutral,mol.e3pt]
            colors='#EA35F8'
        if mol.prefix == 'tri':  
            curves = [mol.neutral,mol.e4pt]
            colors='#367E21'
    
        if region == 'high': xlim = [3500,3200]
        if region == 'low':  xlim = [1700,1300]
         
        curves = [curves] if isinstance(curves, IRSpec) else curves
        if not xlim:
            xlim = [max(mol.freq),min(mol.freq)]
        ax[row,col].tick_params(width=2)
        for axis in ['top','bottom','left','right']:
            ax[row,col].spines[axis].set_linewidth(2)
    
        if interpolate:
            start = curves[0]
            finish = curves[-1]
            #ax.set_prop_cycle('color',plt.cm.Reds(np.linspace(0.3,1,num)))
            for step in np.linspace(0,num,num=num):
                curve = ((num - step)/num)*start.curve + (step/num)*finish.curve
                if step == 0:
                    label = start.name
                    ls = '-'
                    lw = 2 
                    color = 'black'
                    zorder = 2
                elif step == np.linspace(0,num,num=num)[-1]:
                    label = finish.name
                    if not colors:
                        color = 'r'
                    else:
                        color = colors
                    zorder = 3
                    ls = '-'
                    lw = 2 
                else:
                    label = None
                    color = colors 
                    zorder = 1
                    ls = '-'
                    lw = 1
                ax[row,col].plot(mol.freq,curve,label=label,color=color,zorder=zorder,lw=lw)
           #plt.ylim([-50,max(max(start.curve),max(finish.curve))])
            ax[row,col].set_ylim([-50,1200])
            #w = 3
            #plt.axvspan(1659+w, 1638-w, facecolor='r', alpha=0.2,zorder=-1,label='imine C=N str')
            #plt.axvspan(1616+w, 1563-w, facecolor='orange', alpha=0.2,zorder=-1,label='aromatic sym C=C str')
            #plt.axvspan(1532+w, 1502-w, facecolor='y', alpha=0.2,zorder=-1,label='interaromatic C-C str')
            #plt.axvspan(1493+w, 1477-w, facecolor='g', alpha=0.2,zorder=-1,label='benzimidazole$_1$ NH bend + phenol C=O str')
            #plt.axvspan(1471+w, 1376-w, facecolor='b', alpha=0.2,zorder=-1,label='t-butyl C-H bend')
            #plt.axvspan(1428+w, 1366-w, facecolor='purple', alpha=0.2,zorder=-1,label='aromatic asym C=C str')
            #plt.axvspan(1355+w, 1334-w, facecolor='pink', alpha=0.2,zorder=-1,label='benzimidazole C=N str')
    
        ax[row,col].legend()
        ax[row,col].set_xlim(xlim)
        ax[row,col].set_xticks(np.arange(min(xlim),max(xlim),50))
        ax[row,col].set_xticks(np.arange(min(xlim),max(xlim),25),minor=True)
       
        ax[row,col].set_yticks([0])
        #plt.grid(color='gray',alpha=0.4)
        ax[row,col].grid(color='gray',alpha=0.4,which='minor',ls='--')
if save:
    if isinstance(save,str):
        #plt.tight_layout()
        plt.savefig(save,dpi=300,bbox_inches='tight')
    else:
        raise ValueError("'save' needs to be your filename!")
if show:
    plt.show()
plt.close()


