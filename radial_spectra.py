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
CENTER = [160,160]
WL_NUMBER = len(cube[:,1,1].data)
binned_spectra = np.zeros((WL_NUMBER, BIN_NUMBER))
im = np.array(cube[2600,:,:].data)
R = np.arange(0, Rmax+ANN_SIZE, step=ANN_SIZE)


for x,wl in enumerate(cube[:,0,0]):
    ima = cube[x,:,:].data
    
    for i, r in enumerate(R):
        bin_flux = np.array([])
        # print("Running aperture number",i)
        if i == 0:
            aperture = CircularAperture(CENTER, ANN_SIZE).to_mask(method="exact")
        else:
            aperture = CircularAnnulus(CENTER, r, r+ANN_SIZE).to_mask(method="exact")

        reg = aperture.cutout(ima)
        print(reg) # This is always none... why?
        reg_avg = np.nanmedian(reg)
        bin_flux = np.append(bin_flux, reg_avg)
    print(bin_flux)
    binned_spectra[x] = bin_flux
    

        

    