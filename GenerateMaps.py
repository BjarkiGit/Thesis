import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from natsort import natsorted
# from Fitter_working import fit
from Fitter_working_branch_ARMod import fit
from mpdaf.obj import Cube
from astropy.io import fits as f
from models import gaussFit



"""
# For J0156
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/"
DATA = PATH+"J0156_DATACUBE_FINAL.fits"
xRange = np.arange(140, 180, step=1) # J0156
yRange = np.arange(130, 170, step=1) # J0156
Z_INIT = 0.2696 # J0156

# For J0004
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0004/"
DATA = PATH+"J0004_DATACUBE_FINAL.fits"
xRange = np.arange(140, 180, step=1) # J0004
yRange = np.arange(140, 180, step=1) # J0004
Z_INIT = 0.2386 # J0004
"""

# For J0139
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0139/"
DATA = PATH+"J0139_DATACUBE_FINAL.fits"
xRange = np.arange(140, 180, step=1) # J0139
yRange = np.arange(140, 180, step=1) # J0139
Z_INIT = 0.2386 # J0004

lines = pd.read_csv("Lines.txt", sep = r"\s+")
cube = Cube(filename = DATA)
var = f.open(DATA)[2].data
# The guess for z is what I think it "should" be, the FWHM guess is just random
z_val = Z_INIT
WIDVAL = 2.0

files = natsorted(glob(PATH+"results/"+"*.csv"))
# There's probably a better way to set up this loop rather than the nested loops I have now
for x in xRange:
    # I feel like I want to pass in a different z when I go back to the other edge
    # for now it's just the one I had at the start
    z_val = Z_INIT
    for y in yRange:
        name = PATH+"results/fitResult"+str(x)+"_"+str(y)+".csv"
        printname = str(x)+","+str(y)
        if name in files:
            print("Fit of",printname,"exists, passing")
            
        else:
            print("Running pixel ["+printname+"]")
            result, df, fitResult = fit(cube, lines, var, x, y, z_val, WIDVAL)
            if result is not None:
                z = result.params["z"]

                print("Fit successfully")
                # The dataframe doesn't like saving things of different lenghts
                # so I set up a column with one value repeated all the way down
                # to save in the dataframe for z, z_err, FWHM, and FWHM_err

                # The try/except was set up to try to prevent it from breaking when z error
                # can't be estimated
                # This never returns z_val as an array. I've checked that
                length = np.ones_like(fitResult["flux"])
                
                try:
                    zs = length*z.value
                    zerrs = length*z.stderr

                    
                except Exception as e:
                    zerrs=length*999
                    print(e," in z calculation")

                wid = result.params["wid"]
                WIDVAL = np.abs(wid.value)
                ws = length*WIDVAL

                try:
                    werrs = length*wid.stderr

                except Exception as e:
                    werrs = length*999
                    print(e, "in FWHM calculation")
            else:
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
            print(name, "saved")

            # I want to save a plot of all the fits, just haven't figured out how yet
            """
            WL = df["wl"].values

            plt.plot(WL, df["flux0"], label="Flux")
            plt.plot(WL, gaussFit(result.params, WL, lines = lines), label="Fit")
            TITLE = "Pix ("+str(x)+","+str(y)+"), z="+str(z_val)
            plt.title(TITLE)
            plt.legend(loc=(1,0))
            plt.savefig(PATH+"results/figs/test/"+str(x)+"_"+str(y)+".pdf")
            """