import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as f
from lmfit import Model, Parameters, Minimizer, minimize, report_fit, fit_report
import pandas as pd
from mpdaf.obj import Cube
from models import *

path = "../J0156/"

cube = Cube(filename = path+"ADP.2021-10-13T20:49:23.745.fits")
spe = cube[:,160,160]
range = spe.get_range()
step = spe.get_step()

wl = np.arange(range[0], range[1]+step, step)
flux = spe.data


df = pd.DataFrame({"flux0":flux, "wl":wl})
# Total number of wavelengths is 3829
wdw = 100
newFlux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True).median()
df["flux"] = newFlux
nan_count = df["flux"].isna().sum()
print(nan_count)
print(len(df["flux"]))
df = df.fillna(0)


# Getting the lines from a file
lines = pd.read_csv("../Lines.txt", sep = r"\s+")

wl_vac = lines["wl"]
lineName = lines["line"]
st = 2 # Initial guess for strength of Halpha
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

result = minimize(gauss, pars, args=(df["wl"].values,), kws={"f":df["flux"].values, "lines":lines}, nan_policy="omit")

df["fit"] = df["flux"].dropna()+result.residual
print(df.head().interpolate(method="linear"))
print(result.params.pretty_print())
print(f"Max: {np.max(np.abs(result.residual))}, mean: {np.mean(result.residual)}, median: {np.median(result.residual)}")
plt.plot(df["wl"], df["flux"])
plt.plot(df["wl"], (df["flux"].dropna()+result.residual).interpolate(method="linear"))
plt.show()
