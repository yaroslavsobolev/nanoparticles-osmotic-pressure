import numpy as np
import matplotlib.pyplot as plt
import skimage
import scipy
from lmfit.models import LorentzianModel

mod = LorentzianModel()

image = skimage.io.imread('data/MEF_80_20_24h_0005__for_dist.jpg')
nm_per_px = 200/890 #nanopeters per pixel in this image

def get_proj_at_angle(angle):
    image2 = scipy.ndimage.rotate(image, angle=angle)
    dims = image2.shape
    wsize = 150
    hsize = 140
    w0 = int(round(dims[0]/2))+80
    h0 = int(round(dims[1]/2))-30
    image3 = image2[w0-wsize:w0+wsize, h0-hsize:h0+hsize, :]
    projection = np.mean(image3, axis=1)[:,0]
    return projection

good_angle = 37.3169
projection = get_proj_at_angle(good_angle)
plt.show()

fig0 = plt.figure(0)
plt.plot(np.linspace(0, projection.shape[0]*nm_per_px, projection.shape[0]), projection)

# Number of samplepoints
N = projection.shape[0]
# sample spacing
T = nm_per_px
y = projection
yf = scipy.fftpack.fft(y)
xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

fig1, ax = plt.subplots(1)
thegraphy = 2.0/N * np.abs(yf[:N//2])
fromj = 9
toj = 17
ax.plot(xf[:50], thegraphy[:50], '-o')
y = thegraphy[fromj:toj]
x = xf[fromj:toj]
pars = mod.guess(y, x=x)
out = mod.fit(y, pars, x=x)
x2 = np.linspace(np.min(x), np.max(x), 300)
y2 = mod.eval(out.params, x=x2)
plt.plot(x2, y2, '--', alpha=1)
print(out.fit_report(min_correl=0.25))

fromj = 18
toj = 25
y = thegraphy[fromj:toj]
x = xf[fromj:toj]
pars = mod.guess(y, x=x)
out = mod.fit(y, pars, x=x)
x2 = np.linspace(np.min(x), np.max(x), 300)
params2 = out.params
# params2['center'] = 0.15390835*1.3
y2 = mod.eval(params2, x=x2)
plt.plot(x2, y2, '--', alpha=1)
print(out.fit_report(min_correl=0.25))

plt.show()