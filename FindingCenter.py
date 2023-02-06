import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd
from mpdaf.obj import Cube
from models import *


path = "../J0156/"

cube = Cube(filename = path+"ADP.2021-10-13T20:49:23.745.fits")
im = np.array(cube[1199,:,:].data)
spe = cube[:,159,152]


plt.imshow(im, cmap='hot', interpolation='nearest')
plt.show()

spe.plot()
plt.show()
