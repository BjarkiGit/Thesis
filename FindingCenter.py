import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
from models import *
from astropy.io import fits as f
path = "Data/J0156/"
dataFile = path+"J0156_DATACUBE_FINAL.fits"

cube = Cube(filename=dataFile)
hdul = f.open(dataFile)

# im = np.array(cube[1199,:,:].data)
im = hdul[1].data[2980]
spe = cube[:,150,160]


plt.imshow(im, cmap='hot', interpolation='nearest')
plt.show()

spe.plot()
plt.show()
