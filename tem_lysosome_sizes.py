import numpy as np
import matplotlib.pyplot as plt

dtypes = ['str']
dtypes.extend([np.dtype(np.float64)]*5)
dtypes = np.array(dtypes)
def load_file(filename):
    return np.loadtxt(filename, delimiter='\t', skiprows=2, usecols=[1,2,3,4,5])
ht = load_file('data/HT8020_24h.txt')
mda = load_file('data/MDA8020_24h.txt')
mef = load_file('data/MEF8020_24h.txt')

def func(x):
    return np.sqrt(x[:,0]/np.pi)*2/1000

fig,ax = plt.subplots(1, figsize=(2.3,4))
ax2 = ax.twinx()
bp = ax.boxplot([func(x) for x in [ht, mda, mef]],
                 0, 'rs', 1, whis='range', patch_artist=True,
                 widths=0.6)
for i, x in enumerate([ht, mda, mef]):
    data = func(x)
    ax.plot(np.random.normal(i+1, 0.1, size=len(data)), data,
         '.', color='green',
         alpha=1, markersize=5, zorder=10)
for i,x in enumerate(bp['boxes']):
    if i == 0 or i == 1:
        x.set_facecolor('magenta')
    else:
        x.set_facecolor('grey')
    x.set_alpha(0.5)
for median in bp['medians']:
    median.set(color='blue', linewidth=3, linestyle='-.')

ax.set_ylabel('Effective diameter of lysosome section, μm')
ax.set_xticklabels(['HT', 'MDA', 'MEF'])
lims = ax.get_ylim()
ax2.set_ylim((4/np.pi*lims[0], 4/np.pi*lims[1]))
ax2.set_ylabel('Corrected diameter of lysosome, μm')
plt.tight_layout()
fig.savefig('tem_effective_diameters.png', dpi=600)
plt.show()