import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from models import *
from Fitter import fit

path = "../J0156/"
dataFile = path+"ADP.2021-10-13T20:49:23.745.fits"

# xRange = np.arange(157, 159, step=1)
# yRange = np.arange(151, 153, step=1)


lines = pd.read_csv("../Lines.txt", sep = r"\s+")

result, df = fit(dataFile, lines, 157, 153, 100, 5, 0.27, 2, 400)


# Plotting the fit
plt.plot(df["wl"], df["maskFlux"], linewidth=0.8)
plt.plot(df["wl"], gaussFit(result.params, df["wl"], lines=lines["wl"]), 
        color = "r", ls="--", linewidth = 0.5)
plt.xlabel(r"Wavelength $\AA$")
plt.ylabel(r"Flux [erg/($\AA$ cm$^2$) x10$^{-20}$]")
plt.show()