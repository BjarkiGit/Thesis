import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpdaf.obj import Cube, iter_spe, iter_ima # iter_spe lets me loop over all spectra, iter_ima lets me loop over all images
from astropy.io import fits as f
from photutils.aperture import aperture_photometry, CircularAnnulus, CircularAperture
from lmfit import Parameters, minimize
from models import gaussFit

names = ("J0004","J0156","J0139","J0232","J2318")

NAME = names[1]

PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/"
DATA = PATH+NAME+"_DATACUBE_FINAL.fits"
cube = Cube(filename = DATA)[:,100:220,100:220]


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

BIN_NUMBER = 10
ANN_SIZE = Rmax/BIN_NUMBER
WL_NUMBER = len(cube[:,1,1].data)
binned_spectra = np.zeros((WL_NUMBER, BIN_NUMBER+1))
# binned_spectra2 = np.zeros((WL_NUMBER, BIN_NUMBER+1))
image = np.array(cube[2600,:,:].data)
CENTER = np.unravel_index(np.nanargmax(image), image.shape)[::-1]
R = np.arange(0, Rmax+ANN_SIZE, step=ANN_SIZE)

spe = cube[:,CENTER[1],CENTER[0]]
wl_range = spe.get_range()
step = spe.get_step()
z_corrected_wl = np.arange(wl_range[0], wl_range[1]+step, step)
""" 
for wl, dummy in enumerate(cube[:,0,0]):
    ima = cube[wl,:,:].data
    bin_flux = np.array([])
    # bin_flux2 = np.array([])

    for r in R:
        if r == 0:
            aperture = CircularAperture(CENTER, ANN_SIZE).to_mask(method="exact")

        else:
            aperture = CircularAnnulus(CENTER, r, r+ANN_SIZE).to_mask(method="exact")

        # print(aperture)
        reg = aperture.cutout(ima)
        # print(reg) # This is always none... why?
        reg_avg = np.nanmedian(reg)
        # reg_avg2 = np.nanmean(reg)
        bin_flux = np.append(bin_flux, reg_avg)
        # bin_flux2 = np.append(bin_flux2, reg_avg2)
        

    binned_spectra[wl] = bin_flux
    # binned_spectra2[wl] = bin_flux2

    

np.save(PATH+"binnedSpectraTest/Binned", binned_spectra)
# np.save(PATH+"binnedSpectraTest/Binned2", binned_spectra2)
 """


binned_spectra = np.load(PATH+"binnedSpectraTest/Binned.npy")
# binned_spectra = np.load(PATH+"binnedSpectraTest/Binned2.npy")
print(np.shape(binned_spectra))


for spectrum in range(binned_spectra.shape[1]):
    lab = str(R[spectrum])
    plt.plot(z_corrected_wl, binned_spectra[:,spectrum], label = lab)

plt.title(NAME+" median radial spectra")
plt.legend(title = "Outer radius of\n bin in pixels")
plt.show()


# plt.plot(wl, binned_spectra[:,7])
# plt.show()


# for spectrum in range(binned_spectra2.shape[1]):
#     lab = "Bin no. "+str(spectrum)
#     plt.plot(binned_spectra2[:,spectrum], label = lab)

# plt.legend()
# plt.show()

