import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd
from mpdaf.obj import Cube
from models import *


def fit(dataFile, linefile, xPix, yPix, medwdw, stHa, zg, wg):
	cube = Cube(filename = dataFile)
	spe = cube[:, xPix, yPix]

	# The wavelengths are not saved as an array, but I have the
	# start, stop and, step size
	range = spe.get_range()
	step = spe.get_step()
	wl = np.arange(range[0], range[1]+step, step)
	# Call this flux0 as I will (attempt) to remove the continuum from this
	flux0 = spe.data

	# Setting up the data frame
	df = pd.DataFrame({"flux0":flux0, "wl":wl})

	# Take the rolling median of the flux in the data frame and storing
	wdw = medwdw # Window size of rolling median
	flux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True).median()
	df["flux"] = flux
	df = df.fillna(0) # Gets rid of na values


	
	lines = linefile
	wl_vac = lines["wl"]
	lineName = lines["line"]

	# Setting up for some masking and initial guesses for fitting
	st = stHa # Initial guess for strength of Halpha
	strength = st * lines["strength"]
	zguess = zg # Inital guess of redshift
	wguess = wg # Initial guess of width of Gaussian

	# A masked flux which starts out with all zeros,
	# the fluxes close to line centers will be added to it to use
	# in fitting
	fluxMask = np.zeros_like(flux)
	for line in wl_vac:
		Center = line*(1+zguess)
		Width = 20*wguess
		""" # This is just for visualisation purposes, not needed if things work
		window = np.arange(Center-Width, Center+Width+step, step)
		plt.fill_between(window, 0, 100, alpha=0.4) """
		# Now I loop over all wavelengths and if they are within "width"
		# of the estimated line center I add them, otherwise the flux is
		# kept at 0
		i = 0
		for wl in df["wl"]:
			if abs(wl-Center) < Width:
				fluxMask[i] = flux[i]
			else:
				pass
			i += 1
	# Saving the masked flux in the dataframe	
	df["maskFlux"] = fluxMask

	# Creating parameters for the fit
	pars = Parameters()

	pars.add_many(("amp0", strength[0], True, 0, 500),
				("amp1", strength[1], True, 0, 500),
				("amp2", strength[2], True, 0, 500),
				("amp3", strength[3], True, 0, 500),
				("amp4", strength[4], True, 0, 500),
				("amp5", strength[5], True, 0, 500),
				("amp6", strength[6], True, 0, 500),
				("amp7", strength[7], True, 0, 500),
				("amp8", strength[8], True, 0, 500),
				("amp9", strength[9], True, 0, 500),
				("amp10", strength[10], True, 0, 500),
				("amp11", strength[11], True, 0, 500),
				("amp12", strength[12], True, 0, 500),
				("amp13", strength[13], True, 0, 500),
				("amp14", strength[14], True, 0, 500),
				("amp15", strength[15], True, 0, 500),
				("amp16", strength[16], True, 0, 500),
				("amp17", strength[17], True, 0, 500))		
	pars.add("z", zguess, True, 0, 1)
	pars.add("wid", wguess, True, 0, 10)

	# Minimizing residual
	result = minimize(gaussFit, pars, args=(df["wl"].values,), kws={"f":df["maskFlux"].values, "lines":wl_vac}, nan_policy="omit")

	# Printing resulting parameters and some statistics on the residual
	""" 	
	print(result.params.pretty_print())
	print(f"Max: {np.max(np.abs(result.residual))}, mean: {np.mean(result.residual)}, median: {np.median(result.residual)}") 
	"""

	return result, df