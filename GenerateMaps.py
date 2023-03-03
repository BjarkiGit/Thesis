import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# from Fitter_working import fit
from Fitter_working_branch_ARMod import fit
from mpdaf.obj import Cube
from models import gaussFit

# PATH = "/Data/J0156/" Don't know why this stopped working...



"""
# For J0156
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/"
DATA = PATH+"J0156_DATACUBE_FINAL.fits"
xRange = np.arange(140, 180, step=1) # J0156
yRange = np.arange(130, 170, step=1) # J0156
Z_INIT = 0.27 # J0156
"""
# For J0004
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0004/"
DATA = PATH+"J0004_DATACUBE_FINAL.fits"
xRange = np.arange(155, 175, step=1) # J0004
yRange = np.arange(150, 170, step=1) # J0004
Z_INIT = 0.24 # J0004


lines = pd.read_csv("Lines.txt", sep = r"\s+")
cube = Cube(filename = DATA)
# The guess for z is what I think it "should" be, the FWHM guess is just random
z_val = Z_INIT
WIDVAL = 2.0

# There's probably a better way to set up this loop rather than the nested loops I have now
for x in xRange:
    # I feel like I want to pass in a different z when I go back to the other edge
    # for now it's just the one I had at the start
    z_val = Z_INIT
    for y in yRange:
        result, df, fitResult = fit(cube, lines, x, y, z_val, WIDVAL)
        z = result.params["z"]

        print(f"[{x},{y}] run successfully")
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
            z_val = z.value
            
        except Exception as e:
            z_val = z.value
            zerrs=length*999
            print(e)

        wid = result.params["wid"]
        WIDVAL = np.abs(wid.value)
        ws = length*WIDVAL

        try:
            werrs = length*wid.stderr
        except:
            werrs = length*999

        # Adding to the fitResult dataframe and saving with the name of the pixel
        fitResult["z"] = zs
        fitResult["zerr"] = zerrs
        fitResult["FWHM"] = ws
        fitResult["FWHMerr"] = werrs
        fitResult.to_csv(PATH+"results/fitResult"+str(x)+"_"+str(y), index=False)

        if fitResult["flux"][0]/fitResult["flux_err"][0] < 0.5:
            z_val = Z_INIT

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
       