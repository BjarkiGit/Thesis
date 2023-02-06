import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd

# This is where the Gaussian lives
from models import *

path = "../Data/SDSS-Spectra/"

hdul = f.open(path+"spec-0679-52177-0011.fits")
data = hdul[1].data


# Setting up my dataframe from the data
wl = 10**(data["loglam"])
flux = data["flux"]
df = pd.DataFrame({"flux0":flux, "wl":wl})
# Total number of wavelengths is 3829
wdw = 100
newFlux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True, min_periods=1).median()
#df = df.dropna(axis=0)
df["flux"] = newFlux



# Getting the lines from a file
lines = pd.read_csv("../Lines.txt", sep = r"\s+")

wl_vac = lines["wl"]
lineName = lines["line"]
st = 100 # Initial guess for strength of Halpha
strength = st * lines["strength"]
line = dict(zip(lineName, wl_vac))

# print(lines.head())	

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
pars.add("z", 0.25, True, 0, 1)
pars.add("wid", 2, True, 0, 10)

# fitter = Minimizer(gauss, pars, fcn_args=(df["wl"].values,),fcn_kws={"f":df["flux"].values, "lines":lines})
# result0 = fitter.minimize(method="leastsq")
# print(result0.params.pretty_print())

result = minimize(gauss, pars, args=(df["wl"].values,), kws={"f":df["flux"].values, "lines":lines})
print(result.params.pretty_print())
print(f"Max: {np.max(result.residual)}, mean: {np.mean(result.residual)}, median: {np.median(result.residual)}")
plt.plot(df["wl"], df["flux"])
plt.plot(df["wl"], df["flux"]+result.residual)
plt.show()

""" 
G = gauss(df["wl"], z = 0.25, wid = 2, amp = strength*st, cen = wl_vac, flux = df["flux"])
pars = Parameters()
pars.add("z", 0.25)
pars.add("amp", st*strength)
pars.add("cen", wl_vac)			
pars.add("wid", 2)
pars["cen"].vary = False
result = minimize(gauss, pars, args=df["wl"], kws=df["flux"])
 """
"""
st = 100
pars = Parameters()
pars.add_many(("z", 0.25),
		#("amp", 200),
		#("cen", 4000),
		("amp", st*strength),
		("cen", wl_vac),			
		("wid", 2))
#ars["cen"].vary = False
print(df["flux"])
result = minimize(gaussTest2, pars, args = (df["wl"]), kws={'data':df["flux"]})
"""
"""
st = 100
mod = Model(gaussTest)
gFit = np.array([])
for i in range(len(strength)):
	a = 100*strength[i]
	c = wl_vac[i]
	params = mod.make_params(z = 0.25,
				amp = a,
				cen = c,
				wid = 2)
	params["cen"].vary = False
	params["amp"].min = 0
	result = mod.fit(df["flux"], params, x = df["wl"])
	#plt.plot(df["wl"],df["flux"], linewidth=0.5, color="gray")
	#plt.plot(df["wl"],result.best_fit, c="r", ls="--", linewidth=0.8)
	#plt.show()

	if len(gFit) == 0:
		gFit = np.append(gFit, result.best_fit)
	else:
		gFit += result.best_fit
plt.plot(df["wl"],df["flux"], linewidth=0.5, color="gray")
plt.plot(df["wl"],gFit, c="r", ls="--", linewidth=0.8)
plt.show()
"""


                                             
"""                                             
# Fitting model and retrieving results                                             
result = mod.fit(df["flux"], pars, x = df["wl"])
resultArray = np.array(list(result.best_values.items()))
print(result.fit_report())


z = result.best_values["z"]
wid = result.best_values["wid"]

# Getting all error values
i = 0
err = np.zeros_like(result.params)
for par in result.params:
	e = result.params[par].stderr
	err[i] = e
	i += 1

# The first error is for z, the second for the sigma
z_err = err[0]
err = np.delete(err, 0)
wid_err = err[0]
err = np.delete(err,0)

# The rest are alternatingly for the amplitude and for the central wavelength (fixed)
f_err = np.array([])
for i in range(len(err)):
	if i % 2 == 0:
		f_err = np.append(f_err, err[i])


# Now we want to calculate the flux and the error on the flux
amps = ("amp0", "amp1", "amp2", "amp3", "amp4", "amp5", "amp6", "amp7", "amp8", "amp9", 
	"amp10", "amp11", "amp12", "amp13", "amp14", "amp15", "amp16", "amp17")
fit_amps = np.array([])
fluxes = ([])
flux_errs = ([])
i = 0
for amp in amps:
	f = result.params[amp]
	fit_amps = np.append(fit_amps, f)
	flux = np.sqrt(2*np.pi)*wid*f
	fluxes = np.append(fluxes, flux)
	flux_err = np.sqrt((f_err[i]/f)**2+(wid_err/wid)**2)*flux
	flux_errs = np.append(flux_errs, flux_err)
	i += 1
	


# Putting parameters from the fit into a dataframe
fit = pd.DataFrame({"line":lines["line"], "wl":wl_vac, "amp":fit_amps, "amp_err":f_err, 
			"flux":fluxes, "flux_err":flux_errs})

fit.to_csv("code_test.csv", index=False)

# Plotting fit and errors
plt.plot(df["wl"],df["flux"], linewidth=0.5, color="gray")
plt.plot(df["wl"],result.best_fit, c="r", ls="--", linewidth=0.8)
plt.errorbar(fit["wl"]*(1+z), fit["amp"], yerr = fit["amp_err"], fmt="none", barsabove = True,
		capsize=1, elinewidth = 0.7, ecolor="k")

amp = result.best_values["amp0"]


# Loop that prints the names of lines next to them so I can see which lines were 
# not being fit properly
for i in range(18):
	plt.text(lines["wl"][i]*(1+z), lines["strength"][i]*amp, lines["line"][i])



plt.xlabel(r"Wavelength [$\AA$]")
plt.ylabel("Flux")

plt.savefig(path+"/Figs/fit.pdf")
plt.show()
"""

""" print(f'parameter names: {mod.param_names}')
print(f'independent variables: {mod.independent_vars}') """