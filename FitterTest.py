import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd
from mpdaf.obj import Cube
from models import *
from Fitter import fit

path = "../J0156/"
dataFile = path+"ADP.2021-10-13T20:49:23.745.fits"

lines = pd.read_csv("../Lines.txt", sep = r"\s+")

result, df = fit(dataFile, lines, 159, 152, 100, 5, 0.27, 2)



# Plotting the fit
plt.plot(df["wl"], df["maskFlux"], linewidth=0.8)
plt.plot(df["wl"], gaussFit(result.params, df["wl"], lines=lines["wl"]), 
            color = "r", ls="--", linewidth = 0.5)
plt.xlabel(r"Wavelength $\AA$")
plt.ylabel(r"Flux [idk yet, it says somewhere]")
plt.show()