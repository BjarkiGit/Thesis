import sys
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
    sys.exit()

BIN_NUMBER = 10
ANN_SIZE = Rmax/BIN_NUMBER
WL_NUMBER = len(cube[:,1,1].data)
binned_spectra = np.zeros((WL_NUMBER+1, BIN_NUMBER+1))
areas = np.array([])

spe = cube[:,50,50] # Just a random spectrum to get the wavelength range
wl_range = spe.get_range()
wlHa = 6562.80*(1+z)
step = spe.get_step()
wls = np.arange(wl_range[0], wl_range[1]+step, step)
difference_array = np.abs(wls-wlHa)
wl_index = difference_array.argmin()
stack_index = np.arange(wl_index-4, wl_index+4)

image = np.mean(cube[stack_index,:,:].data, axis=0)
plt.imshow(image)
CENTER = np.unravel_index(np.nanargmax(image), image.shape)[::-1] # Opposite indexing
plt.scatter(CENTER[0],CENTER[1])
plt.show()
R = np.arange(0, Rmax+ANN_SIZE, step=ANN_SIZE)
print(f"Creating mean annular spectra with Rmax {Rmax} over {BIN_NUMBER} bins\n"+
    f"The annulus size is {ANN_SIZE} pixels")


for wl, dummy in enumerate(cube[:,0,0]):
    if wl % 100 == 0:
        print(np.round(wl/len(wls),3)*100,"% done")
    ima = cube[wl,:,:].data
    bin_flux = np.array([])

    for r in R:
        if r == 0:
            aperture = CircularAperture(CENTER, ANN_SIZE)
            aperturemask = aperture.to_mask(method="exact")

        else:
            aperture = CircularAnnulus(CENTER, r, r+ANN_SIZE)
            aperturemask = aperture.to_mask(method="exact")
            # Saving the area for flux calculation
        # print(aperture)
        reg = aperturemask.cutout(ima)

        try:
            reg_avg = np.nanmean(reg)/aperture.area

        except ValueError:
            reg_avg = np.nan

            # print(ValueError)

        bin_flux = np.append(bin_flux, reg_avg)
        
    binned_spectra[wl] = bin_flux


cols = ["Bin"+str(i) for i in range(BIN_NUMBER+1)]
df = pd.DataFrame(binned_spectra, columns=cols)
wls = np.insert(wls, 0, 0)
df["wl"] = wls


df.to_csv(PATH+"results/annulus/annulus_spectra.csv", index=False)


# binned_spectra = np.load(PATH+"binnedSpectraTest/Binned.npy")



# print(np.shape(binned_spectra))


for spectrum in range(binned_spectra.shape[1]):
    lab = str(R[spectrum])
    plt.plot(wls, binned_spectra[:,spectrum], label = lab)

plt.title(NAME+" median radial spectra")
plt.legend(title = "Outer radius of\n bin in pixels")
plt.show()


