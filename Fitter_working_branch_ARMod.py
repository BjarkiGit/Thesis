import numpy as np
import matplotlib.pyplot as plt
from lmfit import Parameters, minimize
import pandas as pd

from models import *
"""
This function takes in the data and the lines I want to fit as files,
the coordinates in the data and, some initial guesses as well as an
optional maximum value of amplitudes. It can also be given a "windows"
argument to plot the windows where flux is not masked.
It then returns the results of the fit, the dataframe of the input values
and a dataframe with the amplitude and flux of each line along with errors
dataFile = MUSE datacube
linefile = .txt file with a list of lines, wavelengths and relative strengths
xPix, yPix = Pixels in the dataFile where lines are to be fit
zg = Initial guess for z, used for masking lines
wg = Initial guess of FWHM, used for masking lines
aMax = The maximum amplitude. Mostly used in troubleshooting
windows = Whether or not to plot the unmasked windows where we look for lines
"""


#def get_ha_wl_guess(wl, flux, zguess, winwidth=100):

def fit(cube, linefile, xPix, yPix, zg, wg, aMax=50000, windows=False):
    # Using mpdaf.obj to get the spectrum at xPix and yPix
    spe = cube[:, yPix, xPix]

    # The wavelengths are not saved as an array in the cube,
    # but we have the start, stop and, step size so we can make one
    wl_range = spe.get_range()
    step = spe.get_step()
    wl = np.arange(wl_range[0], wl_range[1]+step, step)
    # Call this flux0 as I will (attempt) to remove the continuum from this
    flux0 = spe.data

    # Setting up the data frame
    df = pd.DataFrame({"flux0":flux0, "wl":wl})

    # Take the rolling median of the flux in the data frame and storing
    wdw = 100 # Window size of rolling median
    flux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True).median()
    df["flux"] = flux
    df = df.fillna(0) # Gets rid of na values

    # Getting my lines into a dataframe and setting up a list
    lines = linefile
    wl_air = lines["wl"]
    lineName = lines["line"] # Not currently used

    # Setting up for some masking and initial guesses for fitting
   
    zguess = zg # Inital guess of redshift
    # if type(zg) != float and type(zg) != np.float64:
    #     zguess = 2.7
    
    # It seems that the width becomes really small
    # It should not be able to be less than 

    #wguess = np.max([float(wg), 2])
    wguess = wg
    if wguess < 1.4:
    	wguess = 1.4
    
    # A masked flux which starts out with all zeros,
    # the fluxes close to line centers will be added to it to use
    # in fitting
    fluxMask = np.zeros_like(flux)
    for u, line in enumerate(wl_air):
        Center = line*(1+zguess)
        Width = 300
        # print('\tCenter: ', Center)
        # print('\tZguess: ', zguess)
        # print('\tLine: ', line)
        # print('\tWidth: ', Width)

        fluxMask[np.where(np.abs(df['wl'].values-Center) <= Width)] = flux.values[np.where(np.abs(df['wl'].values-Center) <= Width)] 


        # THis is a long loop so i replaced it with the line above
        """
        for i, wl in enumerate(df["wl"]):
            try:
                if abs(wl-Center) < Width:
                    pass
                    #fluxMask[i] = flux[i]
                    
                else:
                    pass
            except Exception as e:
                plt.figure()
                plt.plot(df['wl'], df["flux0"])

                plt.plot(df['wl'], df["flux"])
                plt.plot(df['wl'], fluxMask)
                #plt.xlim(wl_air[0]-2*wg, wl_air+2*wg)
                plt.show()


                raise(Exception)
        """

        if line == wl_air[0]:
            st = np.max(fluxMask)
            wlHa = df["wl"][np.where(fluxMask == st)[0]].values

            halfMax = np.argmin(np.abs(fluxMask-st/2))
            wg = np.abs(wlHa-df["wl"].values[halfMax])*2
            zguess = (wlHa/line)-1

    # Setting the strength of the lines based on the amplitude of H\alpha
    strength = st * lines["strength"]

    # Saving the masked flux in the dataframe
    df["maskFlux"] = fluxMask
    df = df.fillna(0)
    # Creating parameters for the fit
    amps = ['amp'+str(i) for i in range(18)]
    pars = Parameters()
    for u, amp in enumerate(amps):
        pars.add(amp, strength[u], True, -1000, aMax)
    
    pars.add("z", zg, True, zg-1e-3, zg+1e-3)
    pars.add("wid", wguess, True, 0, 10)
    # Minimizing residual and retrieving results
    print("Fitting")
    result = minimize(gaussFit, pars, args=(df["wl"].values,), kws={"f":df["maskFlux"].values, "lines":wl_air}, nan_policy="omit")
    z = result.params["z"]
    wid = result.params["wid"]
    wid_err = wid.stderr

    if wid.stderr is None:
        wid_err = 1e10


    # print(wid.stderr)
    # Putting amplitudes (and errors) in array, calculating fluxes and, 
    # calculating errors on fluxes and putting in the arrays
    
    amp = np.array([])
    amp_err = np.array([])
    flux = ([])
    flux_err = ([])

    for i, am in enumerate(amps):
        ap = result.params[am]
        a = ap.value
        # This try/except checks if the error was estimated and set's the
        # error to 1e10 if it wasn't
        if ap.stderr is None:
            a_err = 1e10
        else:
            a_err = ap.stderr
        # Checking if the amplitude is negative
        # if so, set amplitude and flux to 0 and errors to 1e10
        if a <= 0:
            a = 0
            a_err = 1e10
            fl = 0
            ferr = 1e10

        else:
            # print(f"a: {a}, a_err: {a_err}, wid: {wid.value}, wid_err: {wid_err}")
            fl = a*wid.value*4*np.sqrt(np.log(2)*np.pi)
            ferr = np.sqrt((a_err/a)**2+(wid_err/wid)**2)*fl

        amp_err = np.append(amp, a_err)
        amp = np.append(amp, a)
        flux = np.append(flux, fl)
        flux_err = np.append(flux_err, ferr)
        # print(type(wid_err))
            # flux[i] = 0
            # print("Flux error could not be estimated\nAssuming no line detected")
        

        i += 1
    # Saving fluxes and amplitudes of each line into a dataframe.
    fitResult = pd.DataFrame({"line":lines["line"], "wl":wl_air, "amp":amp, "amp_err":amp_err, 
			"flux":flux, "flux_err":flux_err})


    # Returning the minimizer.result the dataframe of the data fed in
    # and the amplitude and flux array
    return result, df, fitResult
