import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_color_codes()

mono = [1.11,0.71,0.41]
di = [1.07,0.70,0.54,0.30]
tri = [1.06,0.72,0.53,0.43,0.21]

avg = [np.mean([1.11,1.07,1.06]),np.mean([0.71,0.70,0.72]),np.mean([0.54,0.53]),np.mean([0.43])]


def f(x,a,b,c):
    x = np.asarray(x)
    y = np.zeros_like(x)
    y = a - b*np.log(x+1)
    #y[-1] += c
    return y

xtri = np.arange(0,5)
xdi = np.arange(0,4)
xmono = np.arange(0,3)
plt.plot(xtri,tri,'ro',label='tri-BIP')
plt.plot(xdi,di,'yo',label='di-BIP')
plt.plot(xmono,mono,'bo',label='mono-BIP')

ydata = avg
xdata = range(0,len(ydata)) 

popt, pcov = curve_fit(f, xdata, ydata)
print(popt)
residuals = ydata - f(xdata, *popt)
ss_res = np.sum(residuals**2)
ss_tot = np.sum((ydata-np.mean(ydata))**2)
r_squared = 1 - (ss_res / ss_tot)

#Note R2 invalid for nonlinear regression
print("R2 = ",r_squared)
plt.xticks([0,1,2,3,4],['E0PT','E1PT','E2PT','E3PT','E4PT'])
plt.xlabel('Number $n$ protons transferred (E$n$PT)')
plt.ylabel('Redox Potential (V)')
plt.legend()
plt.grid(alpha=0.4)

#popt = (1.06459574, 0.47305145, 0.13)
print(*popt)
plt.plot(xtri,f(xtri,*popt),'r--')
plt.plot(xdi,f(xdi,*popt),'y--')
plt.plot(xmono,f(xmono,*popt),'b--')
#plt.text(0.75,1.0,r'1.06 - 0.47ln(n+1) - 0.13$\delta_{n,N}$')
#plt.text(0.75,0.95,'(N is maximal protons transfered)')
#plt.savefig('EnPT-fullfit.pdf',dpi=300)
plt.show()


