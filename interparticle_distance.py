import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('particles_distance.csv', skiprows=1, delimiter=',')
print(1)
d = data[:, [-1, -2]]
from sklearn.metrics import pairwise_distances
dists = pairwise_distances(d)
print(1)
dists = dists.flatten()
dists = dists[dists > 0]
plt.hist(dists, bins=200)
plt.show()
dists = dists[dists<10]
print(np.mean(dists))
print(np.std(dists))