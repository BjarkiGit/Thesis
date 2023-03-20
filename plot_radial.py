import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LogNorm
from astropy import constants as const
from glob import glob
from natsort import natsorted
from astropy.cosmology import LambdaCDM

# Defining constants
c = const.c.value/1e3 # speed of light in km/s
H0 = 70 # km/s/Mpc
PX_SCALE = 0.2 # arcsec / px
cosmo = LambdaCDM(H0 = H0, Om0 = 0.3, Ode0 = 0.7)

def kpc_convert(R,z):
    sec_rad = np.pi/(10800*60)
    d_A = cosmo.angular_diameter_distance(z).value*1e3
    return sec_rad*R*PX_SCALE*d_A

names = ("J0004","J0156","J0139","J0232","J2318")

NAME = names[1]




PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/results/annulus/"
files = natsorted(glob(PATH+"bin"+"*.csv"))


Ha = np.array([])
Oiii = np.array([])
O32 = np.array([])
kin = np.array([])
sig = np.array([])
Hb = np.array([])
Oiii_4363 = np.array([])

for fil in files:
    df = pd.read_csv(fil)

    H1 = df.query("line=='HI_6563'")["flux"].values
    O3 = df.query("line=='OIII_5007'")["flux"].values
    O2 = df.query("line=='OII_3726'")["flux"].values + df.query("line=='OII_3729'")["flux"].values
    OIII_4363 = df.query("line=='OIII_4363'")["flux"].values
    Hbeta = df.query("line == 'HI_4861'")["flux"].values

    if O2 == 0:
        rat = np.nan
    else:
        rat = O3/O2

    z = df["z"][0]
    wid = df["FWHM"]

    Ha = np.append(Ha, H1)
    Oiii = np.append(Oiii, O3)
    O32 = np.append(O32, rat)
    kin = np.append(kin, z)
    sig = np.append(sig, wid)
    Hb = np.append(Hb, Hbeta)
    Oiii_4363 = np.append(Oiii_4363, OIII_4363)
bins = np.arange(0,len(files))

r = kpc_convert((np.arange(2,24, step=2)),z)

print(len(r))
print(len(Ha))
plt.clf()
plt.title(NAME)
plt.scatter(r, Ha, label = r"H$\alpha$")
plt.scatter(r, Oiii, label = "O[III]")
plt.yscale("log")
plt.xlabel("Radius [kpc]")
ylab = "Surface brightness\n"+r"[erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2} \times 10^{-18}$]"
plt.ylabel(ylab)
plt.legend()
plt.savefig(PATH+"/maps/HaOiii.pdf")


plt.clf()
plt.title(NAME)
plt.scatter(bins, Ha, label = r"H$\alpha$")
plt.scatter(bins, Oiii, label = "O[III]")
plt.yscale("log")
plt.xlabel("Bin number")
plt.ylabel(ylab)
plt.legend()
plt.savefig(PATH+"/maps/HaOiiiBin.pdf")
 
plt.clf()
plt.title(NAME)
plt.scatter(r, O32, label = "O32")
plt.legend()
plt.savefig(PATH+"/maps/O32.pdf")

plt.clf()
plt.title(NAME)
plt.scatter(r, kin, label = "redshift")
plt.legend()
plt.savefig(PATH+"/maps/redshift.pdf")

# plt.clf()
# plt.title(NAME)
# plt.scatter(bins, sig, label = "FWHM")
# plt.legend()
# plt.savefig(PATH+"/maps/.pdf")

plt.clf()
plt.title(NAME)
plt.scatter(bins, Ha/Hb, label = "Ha/Hb")
plt.legend()
plt.savefig(PATH+"/maps/Ha-Hb.pdf")
