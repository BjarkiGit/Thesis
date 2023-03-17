"""
Similar function as line_fit.py except this one fits for annular spectra 
so there's not as many factors to consider once the data is input
"""
import numpy as np
from lmfit import Parameters, minimize
import pandas as pd

from models import gaussFit


def bad_fit(size):
    output = np.zeros_like(size)
    result = None
    amp = output
    amp_err = output + 1e10
    flux = output
    flux_err = output + 1e10
    return result, amp, amp_err, flux, flux_err


def annulus_fit(wl, flux0, linefile, zg, wg):

    area = flux0[0]
    flux0 = np.delete(flux0, 0)
    wl = wl[1:]
    df = pd.DataFrame({"flux0":flux0, "wl":wl})
    noisy = False


    # Take the rolling median of the flux in the data frame and storing
    wdw = 100 # Window size of rolling median
    flux = df["flux0"]-df["flux0"].rolling(window=wdw, center=True).median()
    df["flux"] = flux
    df = df.fillna(0) # Gets rid of na values

    # Getting my lines into a dataframe and setting up a list
    lines = linefile
    wl_air = lines["wl"]
    lineName = lines["line"]

    zguess = zg
    wguess = max(wg, 1.4)

    fluxMask = np.zeros_like(flux)
    for u, line in enumerate(wl_air):
        Center = line*(1+zguess)

        # Finding the index of the wavlength of Ha
        
        WIDTH = 200
        fluxMask[np.where(np.abs(df['wl'].values-Center) <= WIDTH)] = flux.values[np.where(np.abs(df['wl'].values-Center) <= WIDTH)] 

        # Uncomment the lines in here when you know what to do with the variance
        # in the annular spectra :)
        if line == wl_air[0]:
            # Checking if the spectrum is noisy around Ha
            # difference_array = np.abs(df["wl"]-Center)
            # wl_index = difference_array.argmin()
            st = np.max(fluxMask)
            # var_sum = np.sqrt(np.sum(var[wl_index-WIDTH:wl_index+WIDTH]))
            # flux_sum = np.sum(flux[wl_index-WIDTH:wl_index+WIDTH])
            # SNR = flux_sum/var_sum
            # if SNR < MASK:
            #     noisy = True

    if not noisy:
        strength = st * lines["strength"]

        # Saving the masked flux in the dataframe
        df["maskFlux"] = fluxMask
        df = df.fillna(0)
        # Creating parameters for the fit
        amps = ['amp'+str(i) for i in range(len(lines))]
        pars = Parameters()
        for u, amp in enumerate(amps):
            pars.add(amp, strength[u], True, -1000, 50000)
        
        pars.add("z", zg, True, zg-1e-3, zg+1e-3) # Using zguess here does not work for some reason
        pars.add("wid", wguess, True, 1.4, 10)
        # Minimizing residual and retrieving results
        # print("Fitting")
        try:
            result = minimize(gaussFit, pars, args=(df["wl"].values,), kws={"f":df["maskFlux"].values, "lines":wl_air}, nan_policy="omit")

            wid = result.params["wid"]
            wid_err = wid.stderr

            if wid.stderr is None:
                wid_err = 1e10

            # Putting amplitudes (and errors) in array, calculating fluxes and,
            # calculating errors on fluxes and putting in the arrays
            
            amp = np.array([])
            amp_err = np.array([])
            flux = np.array([])
            flux_err = np.array([])

            for am in amps:
                ap = result.params[am]
                a = ap.value
                # This checks if the error was estimated and set's the
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
                    fl = a*wid.value*4*np.sqrt(np.log(2)*np.pi)/area
                    ferr = np.sqrt((a_err/a)**2+(wid_err/wid)**2)*fl

                amp_err = np.append(amp, a_err)
                amp = np.append(amp, a)
                flux = np.append(flux, fl)
                flux_err = np.append(flux_err, ferr)

        except Exception as e:
            result, amp, amp_err, flux, flux_err = bad_fit(wl_air)

    else:
        result, amp, amp_err, flux, flux_err = bad_fit(wl_air)
        
    # Saving fluxes and amplitudes of each line into a dataframe.
    fitResult = pd.DataFrame({"line":lineName, "wl":wl_air, "amp":amp, "amp_err":amp_err,
			"flux":flux, "flux_err":flux_err})


    # Returning the minimizer.result the dataframe of the input data
    # and the amplitude and flux array
    return result, df, fitResult
