import re
from glob import glob
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
from astropy import constants as const

def prep(array, shape):
    array[array == 0] = np.nan
    return np.rot90(array.reshape(shape))

# PATH = "/Data/J0156/results"
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/results/"
c = const.c.value/1e3 # speed of light in km/s

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
filelist = natsorted(glob(PATH+"fitResult*"))


for i, file in enumerate(filelist):
    df = pd.read_csv(file)
    file = file.replace(PATH,'')
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
z_mask_norm = z_mask - brightest_z
z_mask_v = z_mask_norm*c

np.save(PATH+"/maps/Ha_mask", Ha_mask)
np.save(PATH+"/maps/z_mask", z_mask)
np.save(PATH+"/maps/z_mask_v", z_mask_v)
np.save(PATH+"/maps/Othree_mask", Othree_mask)
np.save(PATH+"/maps/O32", O32)
np.save(PATH+"/maps/OIII_4363_arr", OIII_4363_arr)


"""
heat_map = cmap = plt.cm.get_cmap("hot").copy()
cmap.set_bad('black',1.)
gre = cmap = plt.cm.get_cmap("plasma").copy()
cmap.set_bad('black',1.)


plt.imshow(Ha_mask, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(Ha_mask), vmax=np.nanmax(Ha_mask)))
cbar = plt.colorbar()
cbar.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"H$\alpha$")
plt.savefig(PATH+"figs/Ha.pdf")
plt.savefig(PATH+"figs/Ha.png")
plt.clf()

plt.imshow(Othree_mask, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(Ha_mask), vmax=np.nanmax(Ha_mask)))
cbar2 = plt.colorbar()
cbar2.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"[OIII]")
plt.savefig(PATH+"figs/OIII.pdf")
plt.savefig(PATH+"figs/OIII.png")
plt.clf()

plt.imshow(OIII_4363_arr, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(OIII_4363_arr), vmax=np.nanmax(OIII_4363_arr)))
cbar3 = plt.colorbar()
cbar3.set_label(r"10$^{-20}$ erg cm$^{-2} s^{-1}$")
plt.title(r"[OIII] 4363")
plt.savefig(PATH+"figs/OIII_4363.pdf")
plt.savefig(PATH+"figs/OIII_4363.png")
plt.clf()

plt.imshow(O32, cmap=gre, interpolation="nearest", norm=LogNorm(vmin=0.1, vmax=2))
cbar4 = plt.colorbar()
# cbar4.set_label(r"")
plt.title(r"[OIII]/[OII]")
plt.savefig(PATH+"figs/O32.pdf")
plt.savefig(PATH+"figs/O32.png")
plt.clf()


cmap = plt.cm.get_cmap("coolwarm").copy()
cmap.set_bad('black',1.)
plt.imshow(z_mask_v, cmap=cmap, interpolation="nearest", vmin = -0.00015*c, vmax = 0.00015*c)
cbar5 = plt.colorbar()
cbar5.set_label(r"km s$^{-1}$")
plt.title(r"$z$")
plt.savefig(PATH+"figs/z.pdf")
plt.savefig(PATH+"figs/z.png")
plt.clf()



plt.imshow(z_mask, cmap=cmap, interpolation="nearest", norm=LogNorm(vmin = 2.695e-1, vmax = 2.703e-1))
cbar5 = plt.colorbar()
# cbar5.set_label(r"km s$^{-1}$")
plt.title(r"$z$")
plt.savefig(PATH+"figs/z_log.pdf")
plt.savefig(PATH+"figs/z_log.png")
plt.clf()



plt.imshow(z_mask, cmap=cmap, interpolation="nearest", vmin = 10, vmax = 300)
cbar5 = plt.colorbar()
cbar5.set_label(r"km s$^{-1}$")
plt.title(r"$z$")
plt.savefig(PATH+"figs/z_arm.pdf")
plt.savefig(PATH+"figs/z_arm.png")
plt.clf()
"""