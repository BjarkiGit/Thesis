import numpy as np
from lmfit import Parameters, minimize
import pandas as pd

from models import gaussFit
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
zg = Initial guess for z, also used for masking lines
wg = Initial guess of FWHM, also used for masking lines
"""

def bad_fit(size):
    output = np.zeros_like(size)
    result = None
    amp = output
    amp_err = output + 1e10
    flux = output
    flux_err = output + 1e10
    return result, amp, amp_err, flux, flux_err

def fit(cube, linefile, Var, xPix, yPix, zg, wg, MASK):
    # Using mpdaf.obj to get the spectrum at xPix and yPix
    spe = cube[:, yPix, xPix]
    var = Var[:, yPix, xPix]
    noisy = False
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


    wguess = wg
    if wguess < 1.4:
        wguess = 1.4
   
    # A masked flux which starts out with all zeros,
    # the fluxes close to line centers will be added to it to use
    # in fitting
    fluxMask = np.zeros_like(flux)
    for u, line in enumerate(wl_air):
        Center = line*(1+zguess)

        # Finding the index of the wavlength of Ha
        
        WIDTH = 200
        fluxMask[np.where(np.abs(df['wl'].values-Center) <= WIDTH)] = flux.values[np.where(np.abs(df['wl'].values-Center) <= WIDTH)] 


        if line == wl_air[0]:
            # Checking if the spectrum is noisy around Ha
            difference_array = np.abs(df["wl"]-Center)
            wl_index = difference_array.argmin()
            st = np.max(fluxMask)
            var_sum = np.sqrt(np.sum(var[wl_index-WIDTH:wl_index+WIDTH]))
            flux_sum = np.sum(flux[wl_index-WIDTH:wl_index+WIDTH])
            SNR = flux_sum/var_sum
            if SNR < MASK:
                noisy = True

    # Setting the strength of the lines based on the amplitude of H\alpha
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
                    fl = a*wid.value*4*np.sqrt(np.log(2)*np.pi)
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
    return result, df, fitResult, SNR
