import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import pandas as pd
from glob import glob
import re
from natsort import natsorted
from astropy import constants as const

def prep(array, shape):
    array[array == 0] = None
    return np.rot90(array.reshape(shape))

# path = "/Data/J0156/results"
path = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/results/"
c = const.c.value/1e3

Ha = ([])
Othree = ([])
zed = ([])
zederr = ([])
z_mask = ([])
Ha_mask = ([])
Othree_mask = ([])
OIII_4363_arr = ([])

O32 = ([])
# Reading all fitResults
filelist = natsorted(glob(path+"fitResult*"))


for i, file in enumerate(filelist):
    df = pd.read_csv(file)
    file = file.replace(path,'')
    # I'm pretty sure there is a better way to find these two
    # just didn't have them at time of writing and this seems to work okay
    H1 = df.query("line=='HI_6563'")["flux"].values
    O3 = df.query("line=='OIII_5007'")["flux"].values
    O2 = df.query("line=='OII_3726'")["flux"].values + df.query("line=='OII_3729'")["flux"].values
    OIII_4363 = df.query("line=='OIII_4363'")["flux"].values
    OIII_4363_err = df.query("line=='OIII_4363'")["flux_err"].values

    if O2 == 0:
        rat = np.nan
    else:
        rat = O3/O2
    H1_err = df.query("line=='HI_6563'")["flux_err"].values
    z = df["z"][0]
    zerr = df["zerr"][0]
    # Finding the pixel coordinates in the file name
    num = re.findall(r'\d+', file)
    xPix = int(num[0])
    yPix = int(num[1])

    # Setting up the pixel range for reshaping of arrays
    if i == 0:
        xMin = xPix
        yMin = yPix


    if i == len(filelist)-1:
        xMax = xPix
        yMax = yPix


    Ha = np.append(Ha, H1)
    Othree = np.append(Othree, O3)
    zed = np.append(zed, z)
    zederr = np.append(zederr, zerr)

    if H1/H1_err > 1:
        z_mask = np.append(z_mask, z)
        Ha_mask = np.append(Ha_mask, H1)
        Othree_mask = np.append(Othree_mask, O3)
        O32 = np.append(O32, rat)

    else:
        z_mask = np.append(z_mask, np.nan)
        Ha_mask = np.append(Ha_mask, np.nan)
        Othree_mask = np.append(Othree_mask, np.nan)
        O32 = np.append(O32, np.nan)

    if OIII_4363/OIII_4363_err > 0.5:
        OIII_4363_arr = np.append(OIII_4363_arr, OIII_4363)
    else:
        OIII_4363_arr = np.append(OIII_4363_arr, np.nan)

x = xMax-xMin+1
y = yMax-yMin+1
shape = (x,y)


# Reshaping arrays
Ha_mask = prep(Ha_mask, shape)
z_mask = prep(z_mask, shape)
Othree_mask = prep(Othree_mask, shape)
O32 = prep(O32, shape)
OIII_4363_arr = prep(OIII_4363_arr, shape)

brightest_pixel = np.unravel_index(np.nanargmax(Ha_mask), Ha_mask.shape)


brightest_z = z_mask[brightest_pixel]
z_mask = z_mask - brightest_z
z_mask = z_mask*c


heat_map = cmap = plt.cm.get_cmap("hot").copy()
cmap.set_bad('black',1.)
gre = cmap = plt.cm.get_cmap("plasma").copy()
cmap.set_bad('black',1.)


plt.imshow(Ha_mask, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(Ha_mask), vmax=np.nanmax(Ha_mask)))
cbar = plt.colorbar()
cbar.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"H$\alpha$")
plt.show()

plt.imshow(Othree_mask, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(Ha_mask), vmax=np.nanmax(Ha_mask)))
cbar2 = plt.colorbar()
cbar2.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"[OIII]")
plt.show()

plt.imshow(OIII_4363_arr, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(OIII_4363_arr), vmax=np.nanmax(OIII_4363_arr)))
cbar3 = plt.colorbar()
cbar3.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"[OIII_4363]")
plt.show()


plt.imshow(O32, cmap=gre, interpolation="nearest", norm=LogNorm(vmin=0.1, vmax=2))
cbar4 = plt.colorbar()
cbar4.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"[OIII]/[OII]")
plt.show()



cmap = plt.cm.get_cmap("coolwarm").copy()
cmap.set_bad('black',1.)
plt.imshow(z_mask, cmap=cmap, interpolation="nearest", vmin = -0.00015*c, vmax = 0.00015*c)
cbar5 = plt.colorbar()
cbar5.set_label(r"km s$^{-1}$")
plt.title(r"$z$")
plt.show()