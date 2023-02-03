import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd

path = "../J0156/"

hdul = f.open(path+"ADP.2021-10-13T20:49:23.745.fits")
data = hdul[1].data
data2 = hdul[2].data

print(np.shape(data2))
plt.plot(data[:,180,180])
# plt.plot(data[2000,180,:])
plt.show()