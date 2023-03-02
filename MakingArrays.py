import re
from glob import glob
from natsort import natsorted
import numpy as np
import pandas as pd
from astropy import constants as const

"""
This reads all the .csv files created by GenerateMaps.py
extracts the relevant values and makes them into numpy arrays
which are saved in order to be plotted in PlottingMaps.py
"""

# A function that replaces all 0 with nan and rotates arrays 90 degrees
def prep(array, shape):
    array[array == 0] = np.nan
    return np.rot90(array.reshape(shape))


PATH = "/home/bjarki/Documents/Thesis/Thesis-1/Data/J0156/results/"
c = const.c.value/1e3 # speed of light in km/s, used for kinematics


# Creating empty arrays to contain the values I want
Ha = ([])
Othree = ([])
zed = ([])
zederr = ([])
z_mask = ([])
Ha_mask = ([])
Othree_mask = ([])
OIII_4363_arr = ([])
O32 = ([])

# Reading all fitResults
filelist = natsorted(glob(PATH+"fitResult*"))

"""
This loop loops over all files and puts the relevant values
into arrays. It also finds the x and y coordinates of the files
and saves the minimum and maximum in order to reshape them later
"""
for i, file in enumerate(filelist):
    df = pd.read_csv(file)
    file = file.replace(PATH,'')

    # Obtaining the values
    H1 = df.query("line=='HI_6563'")["flux"].values
    O3 = df.query("line=='OIII_5007'")["flux"].values
    O2 = df.query("line=='OII_3726'")["flux"].values + df.query("line=='OII_3729'")["flux"].values
    OIII_4363 = df.query("line=='OIII_4363'")["flux"].values
    OIII_4363_err = df.query("line=='OIII_4363'")["flux_err"].values

    if O2 == 0:
        rat = np.nan
    else:
        rat = O3/O2
    
    H1_err = df.query("line=='HI_6563'")["flux_err"].values
    z = df["z"][0]
    zerr = df["zerr"][0]
    # Finding the pixel coordinates in the file name
    num = re.findall(r'\d+', file)
    xPix = int(num[0])
    yPix = int(num[1])

    # Setting up the pixel range for reshaping of arrays
    if i == 0:
        xMin = xPix
        yMin = yPix


    if i == len(filelist)-1:
        xMax = xPix
        yMax = yPix


    Ha = np.append(Ha, H1)
    Othree = np.append(Othree, O3)
    zed = np.append(zed, z)
    zederr = np.append(zederr, zerr)

    # Making a mask from the SNR of Halpha
    if H1/H1_err > 1:
        z_mask = np.append(z_mask, z)
        Ha_mask = np.append(Ha_mask, H1)
        Othree_mask = np.append(Othree_mask, O3)
        O32 = np.append(O32, rat)

    # Putting nans for masked values
    else:
        z_mask = np.append(z_mask, np.nan)
        Ha_mask = np.append(Ha_mask, np.nan)
        Othree_mask = np.append(Othree_mask, np.nan)
        O32 = np.append(O32, np.nan)

    # A seperate 0.5 SNR mask for OIII 4363, which is mostly a
    # curiosity for now
    if OIII_4363/OIII_4363_err > 0.5:
        OIII_4363_arr = np.append(OIII_4363_arr, OIII_4363)
    else:
        OIII_4363_arr = np.append(OIII_4363_arr, np.nan)

# Creating the shape using the maxes and mins from the loop
x = xMax-xMin+1
y = yMax-yMin+1
shape = (x,y)


# Reshaping arrays
Ha_mask = prep(Ha_mask, shape)
z_mask = prep(z_mask, shape)
Othree_mask = prep(Othree_mask, shape)
O32 = prep(O32, shape)
OIII_4363_arr = prep(OIII_4363_arr, shape)

# Finding the brightest pixel of H_alpha which we use to make the kinematic z map
brightest_pixel = np.unravel_index(np.nanargmax(Ha_mask), Ha_mask.shape)


brightest_z = z_mask[brightest_pixel]
z_mask_norm = z_mask - brightest_z
z_mask_v = z_mask_norm*c


# Saving all the arrays
np.save(PATH+"/maps/Ha_mask", Ha_mask)
np.save(PATH+"/maps/z_mask", z_mask)
np.save(PATH+"/maps/z_mask_v", z_mask_v)
np.save(PATH+"/maps/Othree_mask", Othree_mask)
np.save(PATH+"/maps/O32", O32)
np.save(PATH+"/maps/OIII_4363_arr", OIII_4363_arr)
