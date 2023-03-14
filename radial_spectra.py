import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from photutils.aperture import aperture_photometry, CircularAnnulus
from lmfit import Parameters, minimize
from models import gaussFit

names = ("J0004","J0156","J0139","J0232","J2318")

NAME = names[1]

PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"_DATACUBE_FINAL.FITS"


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

R = np.arange(0.01, Rmax+ANN_SIZE, step=ANN_SIZE)

for r in R:
    mask = CircularAnnulus(CENTER, r, r+ANN_SIZE).to_mask
    # Idk man