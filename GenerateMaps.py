import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# from Fitter_working import fit
from Fitter_working_branch_ARMod import fit
from mpdaf.obj import Cube
from models import gaussFit

# PATH = "/Data/J0156/" Don't know why this stopped working...
PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/"
DATA = PATH+"J0156_DATACUBE_FINAL.fits"
cube = Cube(filename = DATA)


# For now this is just a semi-random range around the center
xRange = np.arange(140, 180, step=1) # 140, 180 is a good range for J0156
yRange = np.arange(130, 170, step=1) # 130, 170 is a good range for J0156



lines = pd.read_csv("Lines.txt", sep = r"\s+")
# The guess for z is what I think it "should" be, the FWHM guess is just random
ZVAL = 0.27
WIDVAL = 2.0

# There's probably a better way to set up this loop rather than the nested loops I have now
for x in xRange:
    # I feel like I want to pass in a different z when I go back to the other edge
    # for now it's just the one I had at the start
    ZVAL = 0.27
    for y in yRange:
        result, df, fitResult = fit(cube, lines, x, y, ZVAL, WIDVAL)
        z = result.params["z"]

        print(f"[{x},{y}] run successfully")
        # The dataframe doesn't like saving things of different lenghts
        # so I set up a column with one value repeated all the way down
        # to save in the dataframe for z, z_err, FWHM, and FWHM_err

        # The try/except was set up to try to prevent it from breaking when z error
        # can't be estimated
        # This never returns ZVAL as an array. I've checked that
        length = np.ones_like(fitResult["flux"])
        
        try:
            zs = length*z.value
            zerrs = length*z.stderr
            ZVAL = z.value
            
        except:
            ZVAL = z.value
            zerrs=length*999

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
            ZVAL = 0.27

        """
        plt.plot(result["wl"],result["flux0"], label="Flux")
        plt.plot(result["wl"],gaussFit(result.params, lines = lines), label="Fit")
        title = "Pix ("+str(xPix)+","+str(yPix)+"), z="+str(ZVAL)
        plt.title(title)
        plt.legend(loc=(1,0))
        plt.savefig(PATH+"results/figs/"+str(xPix)+"_"+str(yPix))
        """