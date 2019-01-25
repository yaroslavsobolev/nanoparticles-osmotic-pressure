import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

a = 0.5*8.3e-9*100
# eta = 0.15

c_s = 0.290/2*1000 # mol/m^3
deb_len = 0.3e-9 / np.sqrt((c_s) / 1000)  # in nanometers
kappa = 1 / deb_len
kappa = 2.6/a
def get_z_sat(eta):
    R_ws = a * (1 / eta) ** (1 / 3)
    def fun(kappa_star):
        gamma_0 = np.sqrt(1-(kappa/kappa_star)**4)
        f_plus = ((kappa_star*R_ws + 1)/(2*kappa_star)) * np.exp(-1*kappa_star*R_ws)
        f_minus = ((kappa_star * R_ws - 1) / (2 * kappa_star)) * np.exp(1 * kappa_star * R_ws)
        return gamma_0*(-1 + f_plus*np.exp(kappa_star*a)/a + f_minus*np.exp(-1*kappa_star*a)/a) - 4

    sol = optimize.root(fun, kappa, method='hybr')
    if not sol.success:
        print('fail')
    print('kappa:{0:.2e}'.format(kappa))
    print('kappa_star:{0:.2e}'.format(sol.x[0]))
    kappa_star = sol.x[0]
    gamma_0 = np.sqrt(1 - (kappa / kappa_star) ** 4)
    f_plus = ((kappa_star * R_ws + 1) / (2 * kappa_star)) * np.exp(-1 * kappa_star * R_ws)
    f_minus = ((kappa_star * R_ws - 1) / (2 * kappa_star)) * np.exp(1 * kappa_star * R_ws)
    deriv = gamma_0 * (f_plus * (np.exp(kappa_star * a) / a )*(kappa_star - 1/a) + \
                       f_minus * (np.exp(-1 * kappa_star * a) / a)*(-1*kappa_star - 1/a))
    z_sat = deriv*a
    return z_sat, kappa_star
xs = np.linspace(0.01, 0.3)
ys = [get_z_sat(x)[0] for x in xs]
plt.plot(xs, ys)
plt.show()
# print(z_sat)
