import re
from glob import glob
from natsort import natsorted
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd
from astropy import constants as const

PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/results/maps/"
c = const.c.value/1e3 # speed of light in km/s


Ha_mask = np.load(PATH+"Ha_mask.npy")
Othree_mask = np.load(PATH+"Othree_mask.npy")
OIII_4363 = np.load(PATH+"OIII_4363_arr.npy")
z_mask = np.load(PATH+"z_mask.npy")
z_mask_v = np.load(PATH+"z_mask_v.npy")
O32 = np.load(PATH+"O32.npy")
print(OIII_4363)



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

plt.imshow(OIII_4363, cmap=heat_map, interpolation="nearest", norm=LogNorm(vmin=np.nanmin(OIII_4363), vmax=np.nanmax(OIII_4363)))
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



plt.imshow(z_mask, cmap="plasma", interpolation="nearest", norm=LogNorm(vmin = 2.695e-1, vmax = 2.703e-1))
cbar5 = plt.colorbar()
# cbar5.set_label(r"km s$^{-1}$")
plt.title(r"$z$")
plt.savefig(PATH+"figs/z_log.pdf")
plt.savefig(PATH+"figs/z_log.png")
plt.clf()
