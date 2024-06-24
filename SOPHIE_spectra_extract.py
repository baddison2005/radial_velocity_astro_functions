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

path2 = "/Users/u8010412/Desktop/SONG_DATA/TOI-1431_HD201033/Spectra/SOPHIE_spectra/"
os.chdir(path2)

# fits_file = "tessTOI445_lc.fits"
filename = "SOPHIE.2020-01-12T18_19_20.599_s1d_A"
fits_file = filename + '.fits'
planet_name = 'TOI-1431'
filename_short = '2020-01-12T18_19_20.599_s1d_A'
file_out = 'SOPHIE_1D_out/' + planet_name +'_1D_SOPHIE_spec_'+ filename_short +'.txt'




data = fits.getdata(fits_file) # Get the data
hdr = fits.getheader(fits_file) # Get the full header

#print(np.shape(data)) # produces: (5, 51, 2048) = (layer, order, npix)
#print(hdr.get('slit')) # returns the slit number

#Compute average wavelength between ThAr taken before and after science exposure
#calculate spectral flux by dividing the extracted spectrum by the blaze function
#Calculate flux error by taking intensity and dividing by sqrt of number of spectral elements
length = hdr['NAXIS1']
wavelength = np.empty([length])
flux = np.empty([length])
flux_error = np.empty([length])

#print(data)
wavelength_start = hdr['CRVAL1']
delta_wavelength = hdr['CDELT1']
#print(hdr['CRVAL1'])

k = 0
for i in range(length):
    wavelength[k] = wavelength_start + (i*delta_wavelength)
    flux[k] = abs(data[i])
    flux_error[k] = abs(data[i])/(math.sqrt(length)*abs(data[i]))
    k = k + 1
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

plt.savefig('SOPHIE_1D_out/' + planet_name + '_1D_SOPHIE_spec_' + filename_short + '.pdf',dpi=300)



#Output spectrum as text file
output_spectrum = open(file_out, 'w')
for i in range(len(wavelength)):
    if np.isnan(flux[i]) or np.isnan(flux_error[i]):
        continue
    else:
        output_spectrum.write('{:<12.8f} {:<8.7f} {:<11.10f}\n'.format(wavelength[i], flux[i], flux_error[i]))
#endfor
output_spectrum.close()