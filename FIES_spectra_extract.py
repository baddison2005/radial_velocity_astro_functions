#Author: Brett Addison

from astropy.io import fits
from scipy.stats import sigmaclip
import matplotlib.pyplot as plt
import numpy as np
import math
import os

os.getcwd()

# Change to working directory. Modify the following to your data path
# path = "/Users/u8010412/Documents/Repositories/GitHub/IDL/EXOFASTv2/data_fits/TOI445/" \
#       "MAST_2019-05-09T0446/TESS/tessTOI445"

path2 = "/Users/u8010412/Desktop/SONG_DATA/TOI-1431_HD201033/Spectra/TOI-1431_FIES_spectra/"
os.chdir(path2)

# fits_file = "tessTOI445_lc.fits"
filename = "F4_TOI-1431_FIDe210077_2020-05-22T04-42-22.920.spec"
fits_file = filename + '.fits'
planet_name = 'TOI-1431'
filename_short = '2020-05-22T04-42-22.920'
file_out = 'FIES_1D_out/' + planet_name +'_1D_FIES_spec_'+ filename_short +'.txt'




data = fits.getdata(fits_file) # Get the data
hdr = fits.getheader(fits_file) # Get the full header

print(np.shape(data)) # produces: (5, 51, 2048) = (layer, order, npix)
#print(hdr.get('slit')) # returns the slit number

#Compute average wavelength between ThAr taken before and after science exposure
#calculate spectral flux by dividing the extracted spectrum by the blaze function
#Calculate flux error by taking intensity and dividing by sqrt of number of spectral elements
wavelength = np.empty([2062*90])
flux_norm = np.empty([2062*90])
flux = np.empty([2062*90])
blaze = np.empty([2062*90])
flux_spec = np.empty([2062*90])
flux_error = np.empty([2062*90])

k = 0
for i in range(90):
    for j in range(2062):
        wavelength[k] = data[1,i,j]
        flux_norm[k] = abs(data[3,i,j])
        flux[k] = abs(data[0,i,j])/abs(data[2,i,j])
        blaze[k] = abs(data[2,i,j])
        flux_spec[k] = abs(data[0,i,j])
        flux_error[k] = abs(data[0,i,j])/(math.sqrt(2062*90)*abs(data[2,i,j]))
        k = k + 1
    #endfor
#endfor


#print(max(flux))
#print(max(flux_norm))
#print(max(flux_spec))
#print(max(blaze))
#print(max(abs(data[3,*,*])))


#Plot the 1-D spectra
fig, ax = plt.subplots()

# Plot the timeseries in black circles.
ax.plot(wavelength, flux, 'ko',
       markersize=1)
#plt.ylim(0, 0.2)

# Let's label the axes and define a title for the figure.
fig.suptitle(planet_name + " 1D spectrum from " + filename_short)
ax.set_ylabel("Extracted spectrum/blaze")
ax.set_xlabel("Wavelength (A)")

plt.savefig('FIES_1D_out/' + planet_name + '_1D_FIES_spec_' + filename_short + '.pdf',dpi=300)



#Output spectrum as text file
output_spectrum = open(file_out, 'w')
for i in range(len(wavelength)):
    if np.isnan(flux[i]) or np.isnan(flux_error[i]):
        continue
    else:
        output_spectrum.write('{:<12.8f} {:<8.7f} {:<11.10f}\n'.format(wavelength[i], flux[i], flux_error[i]))
#endfor
output_spectrum.close()