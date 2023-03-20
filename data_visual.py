"""
This code doesn't really do much, just used to troubleshoot
spectra and region selection
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpdaf.obj import Cube, Image
from astropy.io import fits as f


names = ("J0004","J0156","J0139","J0232","J2318")

NAME = names[0]

PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/results/annulus/"


PATH2 = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/"
DATA = PATH2+NAME+"_DATACUBE_FINAL.fits"

cube = Cube(filename = DATA)[:,100:220,100:220]

image = np.array(cube[3000,:,:].data)
CENTER = np.unravel_index(np.nanargmax(image), image.shape)[::-1]

# PATH = "Data/J0232/"
# DATA = PATH+"J0232_DATACUBE_FINAL.fits"

# cube = Cube(filename=DATA)
# hdul = f.open(DATA)

# var = hdul[2].data
# spe = cube[:,150,175]
# val = hdul[1].data[:,150,175]


# wl_range = spe.get_range()
# step = spe.get_step()
# wl = np.arange(wl_range[0], wl_range[1]+step, step)
# # im = np.array(cube[1199,:,:].data)

# im = hdul[1].data[2700]

# """ 
# plt.plot(wl, spe.data, label = "mpdaf", linewidth=0.5, color="k")
# plt.plot(wl, val, label = "astropy", ls="--", linewidth=0.5, color="red")
# plt.vlines(3726*(1+0.2386), ymin=-10, ymax=20, label = "OII 3726")
# plt.legend()
# plt.xlim(3700*(1+0.2386), 5200)
# plt.ylim(-10, 20)

# plt.show()
# """
# # print(f"Median of astropy: {np.nanmedian(val)} and mpdaf: {np.nanmedian(spe.data)}")
# # print(f"Mean of astropy: {np.nanmean(val)} and mpdaf: {np.nanmean(spe.data)}")
# # print(f"Len astropy: {len(val)}, len mpdaf: {len(spe.data)}")
# # print(f"nan astropy: {np.count_nonzero(np.isnan(val))}",
# #       f"nan mpdaf: {np.count_nonzero(np.isnan(spe.data))}")

# # print(np.sum(val==spe.data))


# plt.imshow(im, cmap='hot', interpolation='nearest')
# plt.show()


spectra = pd.read_csv(PATH+"annulus_spectra")

bin0 = spectra["Bin2"].values
wl = spectra["wl"].values

bin0 = np.delete(bin0, 0)
wl = np.delete(wl, 0)

df = pd.DataFrame({"flux0":bin0, "wl":wl})
df = df.fillna(0)
wdw = 100 # Window size of rolling median
flux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True).median()
df["flux"] = flux

plt.plot(df["wl"], df["flux"])
plt.show()

