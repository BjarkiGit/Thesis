"""
Radial binning
"""
import numpy as np
import matplotlib.pyplot as plt
from photutils.aperture import aperture_photometry, CircularAnnulus
from astropy import constants as const
from astropy.cosmology import LambdaCDM

# Defining constants
c = const.c.value/1e3 # speed of light in km/s
H0 = 70 # km/s/Mpc
PX_SCALE = 0.2 # arcsec / px
cosmo = LambdaCDM(H0 = H0, Om0 = 0.3, Ode0 = 0.7)
# NAME = "J0004"
# NAME = "J0156"
# NAME = "J0139"
# NAME = "J0232"
# NAME = "J2318"

names = ("J0004","J0156","J0139","J0232","J2318")

def Area(r_in, r_out):
    return np.pi*(r_out**2-r_in**2)

def kpc_convert(R,z):
    sec_rad = np.pi/(10800*60)
    d_A = cosmo.angular_diameter_distance(z).value*1e3
    return sec_rad*R*PX_SCALE*d_A

def flux_binning(r, map):
    annulus = CircularAnnulus(CENTER, r, r+ANN_SIZE)
    total_flux = aperture_photometry(map, annulus)["aperture_sum"]
    avg_flux = total_flux/Area(r,r+ANN_SIZE)
    return avg_flux*PX_SCALE**2*1e2


for NAME in names:

    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/results/maps/"


    if NAME  == "J0004":
        z = 0.2386 # J0004
        Rmax = 20

    elif NAME == "J0139":
        z = 0.3073 # J0139
        Rmax = 20

    elif NAME == "J0156":
        z = 0.2696 # J0156
        Rmax = 20

    elif NAME == "J0232":
        z = 0.3095 # J0232
        Rmax = 20

    elif NAME == "J2318":
        z = 0.2517 # J2318
        Rmax = 40

    else:
        print("No data for galaxy "+NAME)
        exit()

    print("For "+NAME+" Rmax = "+str(kpc_convert(Rmax, z)))
    BIN_NUMBER = 10
    ANN_SIZE = Rmax/BIN_NUMBER
    Ha_map = np.load(PATH+"Ha_mask.npy")
    Ha_map[np.isnan(Ha_map)] = 0

    OIII_map = np.load(PATH+"Othree_mask.npy")
    OIII_map[np.isnan(OIII_map)] = 0

    O32_map = np.load(PATH+"O32.npy")
    O32_map[np.isnan(O32_map)] = 0

    CENTER = np.unravel_index(np.nanargmax(Ha_map), Ha_map.shape)[::-1]
    R = np.arange(0.1, Rmax+ANN_SIZE, step=ANN_SIZE)

    Ha_binned_flux = np.array([])
    OIII_binned_flux = np.array([])
    O32_binned = np.array([])
    for r in R:
        Ha_binned_flux = np.append(Ha_binned_flux, flux_binning(r, Ha_map))
        OIII_binned_flux = np.append(OIII_binned_flux, flux_binning(r, OIII_map))
        O32_binned = np.append(O32_binned, flux_binning(r, O32_map))



    plt.clf()
    plt.scatter(kpc_convert(R+ANN_SIZE/2, z), Ha_binned_flux, label = r"H$\alpha$")
    plt.scatter(kpc_convert(R+ANN_SIZE/2, z), OIII_binned_flux, label = r"O[III]")
    plt.xlabel("Radius [kpc]")
    plt.ylabel(r"Flux [erg s$^{-1}$ cm$^{-2}$ arcsec$^{-2} \times 10^{-18}$]")
    plt.legend()
    plt.title("Surface brightness of spectral lines for "+NAME)
    plt.yscale("log")
    plt.savefig("/home/bjarki/Documents/Thesis/Thesis-1/Data/Profiles/"+NAME+".pdf")
    plt.savefig(PATH+"figs/radial_profile.pdf")
    plt.savefig(PATH+"figs/radial_profile.png")

    plt.clf()
    plt.xlabel("O32 ratio")
    plt.ylabel("Radius [kpc]")
    plt.title("O32 ratio for "+NAME)
    plt.scatter(kpc_convert(R+ANN_SIZE/2, z), O32_binned)
    plt.savefig("/home/bjarki/Documents/Thesis/Thesis-1/Data/Profiles/"+NAME+"_O32.pdf")
    plt.savefig(PATH+"figs/O32_radial_profile.pdf")
    plt.savefig(PATH+"figs/O32_radial_profile.png")