import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from natsort import natsorted
# from Fitter_working import fit
from annular_fit import annulus_fit
from mpdaf.obj import Cube
from astropy.io import fits as f
from models import gaussFit
from loguru import logger


names = ("J0004","J0156","J0139","J0232","J2318")

NAME = names[1]


if NAME  == "J0004":
    Z_INIT = 0.2386 # J0004
    MASK = 1.5

elif NAME == "J0139":
    Z_INIT = 0.3073 # J0139
    MASK = 1.5

elif NAME == "J0156":
    Z_INIT = 0.2696 # J0156
    MASK = 1

elif NAME == "J0232":
    Z_INIT = 0.3095 # J0232
    MASK = 1

elif NAME == "J2318":
    Z_INIT = 0.2517 # J2318
    MASK = 1

else:
    print("No data for galaxy "+NAME)
    sys.exit()



logger.add("logs/"+NAME+"_annulus_{time}.log")
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/"+NAME+"/results/annulus/"

spectra = pd.read_csv(PATH+"annulus_spectra")

wl = spectra["wl"].values
spectra = spectra.drop(columns=["wl"])

lines = pd.read_csv("Lines.txt", sep = r"\s+")

Z_VAL = Z_INIT
WIDVAL = 2.0




files = natsorted(glob(PATH+"bin"+"*.csv"))
if len(files) == 0:
    print("No spectra found")
    sys.exit()

for i, bins in enumerate(spectra):
    name = PATH+"fitResult_bin_"+str(i)+".csv"
    if name in files:
        print(print("Fit of bin",str(i),"exists, passing"))
    else:
        print("Fitting bin",str(i))
        flux0 = spectra[bins].values
        result, df, fitResult = annulus_fit(wl, flux0, lines, Z_VAL, WIDVAL)

        if result is not None:
                z = result.params["z"]

                print("Fit successfully")
                # The dataframe doesn't like saving things of different lengths
                # so I set up a column with one value repeated all the way down
                # to save in the dataframe for z, z_err, FWHM, and FWHM_err

                length = np.ones_like(fitResult["flux"])
                
                try:
                    zs = length*z.value
                    zerrs = length*z.stderr
 
                except Exception as e:
                    zerrs=length*999
                    logger.info("bad bin",str(i))


                wid = result.params["wid"]
                ws = length*wid

                try:
                    werrs = length*wid.stderr

                except Exception as e:
                    werrs = length*999

        else:
            logger.info("Bad fit for bin",str(i))
            length = np.ones_like(fitResult["flux"])
            zerrs = length*999
            werrs = length*999
            zs = length*0
            ws = length*0

            # Adding to the fitResult dataframe and saving with the name of the pixel
        fitResult["z"] = zs
        fitResult["zerr"] = zerrs
        fitResult["FWHM"] = ws
        fitResult["FWHMerr"] = werrs
        fitResult.to_csv(PATH+"bin"+str(i)+".csv", index=False)