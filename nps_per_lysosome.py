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
    volfrac_of_close_packing = 0.74048
    corrected_radius = 4 / np.pi * np.sqrt(x[:, 0] / np.pi)
    corrected_volume = 4/3*np.pi*corrected_radius**3
    volfrac = x[:, -1]
    np_rad = 8.3/2 # in nm
    np_volume = 4/3*np.pi*np_rad**3
    return corrected_volume*volfrac/np_volume/1e6

np.savetxt('nps_per_lysosome_ht.txt', func(ht))
np.savetxt('nps_per_lysosome_mda.txt', func(mda))
np.savetxt('nps_per_lysosome_mef.txt', func(mef))

fig,ax = plt.subplots(1, figsize=(2,4))
bp = plt.boxplot([func(x) for x in [ht, mda, mef]],
                 0, 'rs', 1, whis='range', patch_artist=True,
                 widths=0.6)
for i, x in enumerate([ht, mda, mef]):
    data = func(x)
    plt.plot(np.random.uniform(low=i+1-0.2, high=i+1+0.2, size=len(data)), data,
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

plt.ylabel('NPs per lysosome, $10^6$')
ax.set_xticklabels(['HT', 'MDA', 'MEF'])
plt.tight_layout()
fig.savefig('np_per_lyso.png', dpi=600)
fig.savefig('np_per_lyso.eps')

fig2,ax = plt.subplots(1, figsize=(1.8,4))
bp = plt.boxplot([func(x) for x in [ht, mef]],
                 0, 'rs', 1, whis='range', patch_artist=True,
                 widths=0.6)
for i, x in enumerate([ht, mef]):
    data = func(x)
    print('{0}: N={1}'.format(i, len(data)))
    plt.plot(np.random.uniform(low=i+1-0.2, high=i+1+0.2, size=len(data)), data,
         '.', color='green',
         alpha=1, markersize=5, zorder=10)
for i,x in enumerate(bp['boxes']):
    if i == 0:
        x.set_facecolor('magenta')
    else:
        x.set_facecolor('grey')
    x.set_alpha(0.5)
for median in bp['medians']:
    median.set(color='blue', linewidth=3, linestyle='-.')

plt.ylabel('NPs per lysosome, $10^6$')
ax.set_xticklabels(['HT', 'MEF'])
plt.tight_layout()
fig2.savefig('np_per_lyso_ht_mef.png', dpi=600)
fig2.savefig('np_per_lyso_ht_mef.eps')
plt.show()