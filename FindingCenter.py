
import matplotlib.pyplot as plt
from mpdaf.obj import Cube

from astropy.io import fits as f
PATH = "Data/J0004/"
DATA = PATH+"J0004_DATACUBE_FINAL.fits"

cube = Cube(filename=DATA)
hdul = f.open(DATA)

# im = np.array(cube[1199,:,:].data)
im = hdul[1].data[2700]
spe = cube[:,165,160]


plt.imshow(im, cmap='hot', interpolation='nearest')
plt.show()

spe.plot()
plt.show()
