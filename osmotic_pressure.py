# -*- coding: utf-8 -*-
"""
Created on Sat Jun 10 00:17:44 2017

@author: Yaroslav I. Sobolev
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

charge_density = 4.7 # ligands per nm^2
k_b = 1.38e-23 # Boltzmann constant, J/K
N_A = 6.02214e23 # Avogadro's number
Bjerr_len = 0.7e-9 # in water, in nanometers
T = 273 + 37    # temperature in Kelvin
c_s = 0.290/2*1000 # mol/m^3
eta = 0.3
a = 0.5*8.3e-9 # with ligands
ligand_length = 1.5e-9
mixfactor = 0.8
Q = mixfactor*4*np.pi*(a-ligand_length)**2*(4.7/(1e-9)**2) # total charge of particle
do_donnan = True

Vnp = 4/3*np.pi*a**3
print(a)
print(a/(eta)**(1/3))
print(0.3e-9/np.sqrt((c_s)/1000)) # in nanometers)
print('Yay!')

def renormed_charge(R_ws, a, kappa):
    def fun(kappa_star):
        gamma_0 = np.sqrt(1-(kappa/kappa_star)**4)
        f_plus = ((kappa_star*R_ws + 1)/(2*kappa_star)) * np.exp(-1*kappa_star*R_ws)
        f_minus = ((kappa_star * R_ws - 1) / (2 * kappa_star)) * np.exp(1 * kappa_star * R_ws)
        return gamma_0*(-1 + f_plus*np.exp(kappa_star*a)/a + f_minus*np.exp(-1*kappa_star*a)/a) - 4

    sol = optimize.root(fun, kappa, method='hybr')
    if not sol.success:
        print('fail')
    kappa_star = sol.x[0]
    gamma_0 = np.sqrt(1 - (kappa / kappa_star) ** 4)
    f_plus = ((kappa_star * R_ws + 1) / (2 * kappa_star)) * np.exp(-1 * kappa_star * R_ws)
    f_minus = ((kappa_star * R_ws - 1) / (2 * kappa_star)) * np.exp(1 * kappa_star * R_ws)
    deriv = gamma_0 * (f_plus * (np.exp(kappa_star * a) / a )*(kappa_star - 1/a) + \
                       f_minus * (np.exp(-1 * kappa_star * a) / a)*(-1*kappa_star - 1/a))
    z_sat = deriv*a*a/Bjerr_len
    return z_sat, kappa_star

def get_pressure_all_counterions(eta, chargeratio):
    Vnp = 4 / 3 * np.pi * a ** 3
    NP_concentration = eta / Vnp / N_A
    c_eff = chargeratio*NP_concentration*(4*np.pi*(a-ligand_length)**2)*(charge_density/(1e-9)**2)
    return k_b * T * c_eff * N_A

def get_pressure(eta, c_s, renorm_osmotic = True):
    # get concentration of NP counterions
    Vnp = 4/3*np.pi*a**3
    NP_concentration = eta/Vnp/N_A
    counterions_concentration = NP_concentration*Q
    if do_donnan:
        counterions_concentration = 0
    
    # get inverse Debye length
    deb_len = 0.3e-9/np.sqrt((c_s + counterions_concentration/2)/1000) # in nanometers
    kappa = 1/deb_len
    
    #get ES cell radius
    Rws = a/(eta)**(1/3)
    #get saturation (renormalized) charge
    if not renorm_osmotic:
        Qsat = 4*a/Bjerr_len*(1+kappa*a)#*(1 + 7.3*eta**2)
        if Qsat > Q:
            Qsat = Q
        # Qsat = a/Bjerr_len*17
        print('Qsat:{0}'.format(Qsat))

        # reduced potential at the edge of Wigner-Seitz cell
        phi_Rws = Qsat*Bjerr_len*np.exp((a-Rws)/deb_len)/Rws/(1 + a*kappa)
        if do_donnan:
            P = 4*k_b*T*(c_s)*N_A*(np.sinh(phi_Rws/2))**2
        else:
            P = 4 * k_b * T * (c_s + counterions_concentration / 2) * N_A * (np.sinh(phi_Rws / 2)) ** 2 + \
                    k_b * T * counterions_concentration * N_A
            print('counterions: {0}'.format(k_b * T * counterions_concentration* N_A))
    else:
        z_eff, kappa_star = renormed_charge(Rws, a, kappa)
        P = k_b * T *(kappa_star**2 - kappa**2)/(4*np.pi*Bjerr_len)
        print('Z_eff:{0}'.format(z_eff))
    return P

xs = np.linspace(0.001, 0.9, num=100)
fig, (ax, ax2) = plt.subplots(2, sharex=True, sharey=False, figsize=(5,5),
                              gridspec_kw=dict(height_ratios=(5,2)))
fig.subplots_adjust(hspace=0.05)
ax.set_ylim([0.0001,10])
ax.set_xlim([0, 0.74048])
ax.set_yscale( "log" )
y1s = np.array([get_pressure(x, 0.320/2*1000) for x in xs])
y2s = np.array([get_pressure(x, 0.280/2*1000) for x in xs])
ax.fill_between(xs, y1s/(101325), y2s/(101325), color = 'b', alpha = 0.5, label='renorm')
ax.set_ylabel('Osmotic pressure across\nthe lysosome membrane, atm')
plt.xlabel('Volume fraction of AuNPs ($\chi_{TMA}/\chi_{MUA}=80:20$)')
ax.axhline(dashes=[3,3], y=1.4, color='grey')
ax.axhline(dashes=[1,3], y=0.003, color='grey')
ht8020_6 = np.loadtxt('HT8020_6h.txt')
ht8020_24 = np.loadtxt('HT8020_24h.txt')
mda8020_24 = np.loadtxt('MDA8020_24h.txt')
medians = [np.median(x) for x in [mda8020_24, ht8020_24, ht8020_6]]
y_medians = np.array([get_pressure(x, 0.300/2*1000) for x in medians])
ax.plot(medians, y_medians/(101325), 'o', color='purple', alpha=0.6)
bp = ax2.boxplot([mda8020_24, ht8020_24, ht8020_6],
                 0, 'rs', 0, whis='range', patch_artist=True,
                 widths=0.6)
ax2.plot(mda8020_24, np.random.normal(1, 0.1, size=len(mda8020_24)),
         '.', color='darkcyan',
         alpha=0.7)
ax2.plot(ht8020_24, np.random.normal(2, 0.1, size=len(ht8020_24)),
         '.', color='darkcyan',
         alpha=0.7)
ax2.plot(ht8020_6, np.random.normal(3, 0.1, size=len(ht8020_6)),
         '.', color='darkcyan',
         alpha=0.7)
for x in bp['boxes']:
    x.set_facecolor('grey')
    x.set_alpha(0.5)
for median in bp['medians']:
    median.set(color='purple', linewidth=3, linestyle='-.')
y0 = -0.2
ax2.plot([0.17, 0.58], [y0, y0], '|', markersize=10, color='darkorange')
ax2.plot([0.17, 0.58], [y0, y0], color='darkorange')
ax2.plot([0.02], [-1], 'D', markersize=6)
ax2.set_ylim([-1.6,3.6])
plt.tight_layout()
fig.savefig('osmotic.png', dpi=300)
fig.savefig('osmotic.eps', dpi=300)
plt.show()