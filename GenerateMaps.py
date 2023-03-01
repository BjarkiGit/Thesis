import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from models import *
# from Fitter_working import fit
from Fitter_working_branch_ARMod import fit
from mpdaf.obj import Cube
from models import gaussFit

# path = "/Data/J0156/" Don't know why this stopped working...
path = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/"
dataFile = path+"J0156_DATACUBE_FINAL.fits"
cube = Cube(filename = dataFile)


# For now this is just a semi-random range around the center
xRange = np.arange(140, 170, step=1)
yRange = np.arange(140, 170, step=1)



lines = pd.read_csv("Lines.txt", sep = r"\s+")
# The guess for z is what I think it "should" be, the FWHM guess is just random
zval = 0.27
widval = 2.0

# There's probably a better way to set up this loop rather than the nested loops I have now
for i in range(len(xRange)):
    # I feel like I want to pass in a different z when I go back to the other edge
    # for now it's just the one I had at the start
    zval = 0.27
    for u in range(len(yRange)):
        xPix = xRange[i]
        yPix = yRange[u]
        print(zval)
        result, df, fitResult = fit(cube, lines, xPix, yPix, zval, widval, windows=True)
        z = result.params["z"]
        

        print(f"[{xPix},{yPix}] run successfully") 
        # The dataframe doesn't like saving things of different lenghts
        # so I just set up a column with one value repeated all the way down
        # to save in the dataframe for z, z_err, FWHM, and FWHM_err

        # The try/except was set up to try to prevent it from breaking when z error
        # can't be estimated
        # This never returns zval as an array. I've checked that
        length = np.ones_like(fitResult["flux"])
        
        try:
            zs = length*z.value
            zerrs = length*z.stderr
            zval = z.value
            
        except:
            zval = z.value
            zerrs=length*999

        wid = result.params["wid"]
        widval = np.abs(wid.value)
        ws = length*widval
        try:
            werrs = length*wid.stderr
        except:
            werrs = length*999

        # Adding to the fitResult dataframe and saving with the name of the pixel
        fitResult["z"] = zs
        fitResult["zerr"] = zerrs
        fitResult["FWHM"] = ws
        fitResult["FWHMerr"] = werrs
        fitResult.to_csv(path+"results/fitResult"+str(xPix)+"_"+str(yPix), index=False)
        """ 
        plt.plot(result["wl"],result["flux0"], label="Flux")
        plt.plot(result["wl"],gaussFit(result.params, lines = lines), label="Fit")
        title = "Pix ("+str(xPix)+","+str(yPix)+"), z="+str(zval)
        plt.title(title)
        plt.legend(loc=(1,0))
        plt.savefig(path+"results/figs/"+str(xPix)+"_"+str(yPix))
        """