import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from natsort import natsorted
# from Fitter_working import fit
from line_fit import fit
from mpdaf.obj import Cube
from astropy.io import fits as f
from models import gaussFit
from loguru import logger


# NAME = "J0156"
# NAME = "J0004"
# NAME = "J0139"
# NAME = "J0232"
NAME = "J2318"

logger.add("logs/"+NAME+"_{time}.log")

if NAME  == "J0004":
# For J0004 # Run with SNR mask 1.5
    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0004/"
    DATA = PATH+"J0004_DATACUBE_FINAL.fits"
    xRange = np.arange(140, 180, step=1) # J0004
    yRange = np.arange(140, 180, step=1) # J0004
    Z_INIT = 0.2386 # J0004
    MASK = 1.5

elif NAME == "J0139":
# For J0139 # Run with SNR mask 1.5
    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0139/"
    DATA = PATH+"J0139_DATACUBE_FINAL.fits"
    xRange = np.arange(140, 180, step=1) # J0139
    yRange = np.arange(140, 180, step=1) # J0139
    Z_INIT = 0.3073 # J0139
    MASK = 1.5

elif NAME == "J0156":
    # For J0156 # Run with SNR mask 1
    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/"
    DATA = PATH+"J0156_DATACUBE_FINAL.fits"
    xRange = np.arange(140, 180, step=1) # J0156
    yRange = np.arange(130, 170, step=1) # J0156
    Z_INIT = 0.2696 # J0156
    MASK = 1

elif NAME == "J0232":
    # For J0232 # Run with SNR mask 1
    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0232/"
    DATA = PATH+"J0232_DATACUBE_FINAL.fits"
    xRange = np.arange(140, 180, step=1) # J0232
    yRange = np.arange(140, 180, step=1) # J0232
    Z_INIT = 0.3095 # J0232
    MASK = 1

elif NAME == "J2318":
    # For J2318 # Run with SNR mask 1
    PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J2318/"
    DATA = PATH+"J2318_DATACUBE_FINAL.fits"
    xRange = np.arange(120, 200, step=1) # J2318
    yRange = np.arange(120, 200, step=1) # J2318
    Z_INIT = 0.2517 # J2318
    MASK = 1

else:
    print("No data for galaxy "+NAME)
    exit()

lines = pd.read_csv("Lines.txt", sep = r"\s+")
cube = Cube(filename = DATA)
var = f.open(DATA)[2].data
# The guess for z is what I think it "should" be, the FWHM guess is just random
Z_VAL = Z_INIT
WIDVAL = 2.0

files = natsorted(glob(PATH+"results/"+"*.csv"))
# There's probably a better way to set up this loop rather than the nested loops I have now
for x in xRange:
    # I feel like I want to pass in a different z when I go back to the other edge
    # for now it's just the one I had at the start
    for y in yRange:
        name = PATH+"results/fitResult"+str(x)+"_"+str(y)+".csv"
        printname = str(x)+","+str(y)
        if name in files:
            print("Fit of",printname,"exists, passing")
            
        else:
            print("Running pixel ["+printname+"]")
            result, df, fitResult, SNR = fit(cube, lines, var, x, y, Z_VAL, WIDVAL, MASK)
            if result is not None:
                z = result.params["z"]


                print("Fit successfully")
                # The dataframe doesn't like saving things of different lenghts
                # so I set up a column with one value repeated all the way down
                # to save in the dataframe for z, z_err, FWHM, and FWHM_err

                # The try/except was set up to try to prevent it from breaking when z error
                # can't be estimated
                # This never returns Z_VAL as an array. I've checked that
                length = np.ones_like(fitResult["flux"])
                
                try:
                    zs = length*z.value
                    zerrs = length*z.stderr
 
                except Exception as e:
                    zerrs=length*999
                    logger.info("Got bad error estimation for pixel: ",
                                 str(x)+","+str(y), "with SNR: ", str(round(SNR,5)))


                wid = result.params["wid"]
                WIDVAL = np.abs(wid.value)
                ws = length*WIDVAL

                try:
                    werrs = length*wid.stderr

                except Exception as e:
                    werrs = length*999

            else:
                logger.info("Bad SNR ("+str(round(SNR,5))+") for pixel "+str(x)+","+str(y))
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
            fitResult.to_csv(name, index=False)
            # print(name, "saved")

            # I want to save a plot of all the fits, just haven't figured out how yet
            """
            WL = df["wl"].values

            plt.plot(WL, df["flux0"], label="Flux")
            plt.plot(WL, gaussFit(result.params, WL, lines = lines), label="Fit")
            TITLE = "Pix ("+str(x)+","+str(y)+"), z="+str(Z_VAL)
            plt.title(TITLE)
            plt.legend(loc=(1,0))
            plt.savefig(PATH+"results/figs/test/"+str(x)+"_"+str(y)+".pdf")
            """