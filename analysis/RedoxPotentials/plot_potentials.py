import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_color_codes()
fig,ax = plt.subplots()

#mono = [1.11,0.71,0.41]
#di = [1.07,0.70,0.54,0.30]
#tri = [1.06,0.72,0.53,0.43,0.21]
tri = np.asarray([122.69,114.90,110.51,108.22,103.23]) - 122.69

xi = np.arange(0,5)
plt.plot(xi,tri,'b-o',label='tri-BIP',markersize=10,lw=3)

plt.xticks([0,1,2,3,4],fontsize=12)#,['zero','one','two','three','four'])
plt.yticks(fontsize=12)#,['zero','one','two','three','four'])
plt.xlabel('Number $n$ protons transferred',fontsize=16)
plt.ylabel('$\Delta G$ (kcal/mol)',fontsize=16)

ax.tick_params(width=2)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(2)
#plt.legend()
plt.grid(alpha=0.4,lw=1.5)
plt.tight_layout()
plt.savefig('gibbs.pdf')
#plt.show()


