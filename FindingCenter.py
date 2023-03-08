"""
This code doesn't really do much, just used to
estimate the location of the galaxy in the frame and 
approximate the redshift from the probable Ha peak
"""
import matplotlib.pyplot as plt
import numpy as np
from mpdaf.obj import Cube, Image
from astropy.io import fits as f


PATH = "Data/J0139/"
DATA = PATH+"J0139_DATACUBE_FINAL.fits"

cube = Cube(filename=DATA)
hdul = f.open(DATA)

var = hdul[2].data
spe = cube[:,150,175]
val = hdul[1].data[:,150,175]


wl_range = spe.get_range()
step = spe.get_step()
wl = np.arange(wl_range[0], wl_range[1]+step, step)
# im = np.array(cube[1199,:,:].data)

im = hdul[1].data[2700]

""" 
plt.plot(wl, spe.data, label = "mpdaf", linewidth=0.5, color="k")
plt.plot(wl, val, label = "astropy", ls="--", linewidth=0.5, color="red")
plt.vlines(3726*(1+0.2386), ymin=-10, ymax=20, label = "OII 3726")
plt.legend()
plt.xlim(3700*(1+0.2386), 5200)
plt.ylim(-10, 20)

plt.show()
"""
# print(f"Median of astropy: {np.nanmedian(val)} and mpdaf: {np.nanmedian(spe.data)}")
# print(f"Mean of astropy: {np.nanmean(val)} and mpdaf: {np.nanmean(spe.data)}")
# print(f"Len astropy: {len(val)}, len mpdaf: {len(spe.data)}")
# print(f"nan astropy: {np.count_nonzero(np.isnan(val))}",
#       f"nan mpdaf: {np.count_nonzero(np.isnan(spe.data))}")

# print(np.sum(val==spe.data))


plt.imshow(im, cmap='hot', interpolation='nearest')
plt.show()
