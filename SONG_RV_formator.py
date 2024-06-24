import pandas as pd
import pyexcel
import numpy as np
import csv
import matplotlib.pyplot as plt
from numpy import mean

directory = '/Users/u8010412/Desktop/SONG_DATA/TOI-1830_HD133725/'
planet_name = 'TOI1830'
input_file = 'HD133725_2020-07-23.txt'
TESS_period = 9.781642
Mid_transit = 1936.906691
sigma_clip = False
sigma_clip_level = 5.0
Exofast_format = True
Radvel_format = True
Exosam_format = True

Mid_transit = Mid_transit + 2457000.0

raw_data = np.genfromtxt(directory + input_file, dtype='str', comments='#', usecols=(2, 4, 6))
data_length = len(raw_data)

Index = np.arange(data_length)
RV_array = pd.DataFrame(index=Index, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])

RV_array['Time'] = raw_data[:,0]
RV_array['RV'] = raw_data[:,1]
RV_array['RV_error'] = raw_data[:,2]
RV_array['Telescope'] = 'SONG'

RV_array['Time'] = RV_array['Time'].astype(float)
RV_array['RV'] = RV_array['RV'].astype(float)
RV_array['RV_error'] = RV_array['RV_error'].astype(float)

#Remove the bad RVs. Sigma clip the bad ones based on the set sigma clip level.

if sigma_clip == True:
    standard_dev = RV_array.RV.std()
    mean_RV = RV_array.RV.mean()
    #remove data below sigma clip level
    RV_array = RV_array.drop(RV_array[RV_array.RV <= mean_RV - (standard_dev*sigma_clip_level)].index)
    #remove data above sigma clip level
    RV_array = RV_array.drop(RV_array[RV_array.RV >= mean_RV + (standard_dev*sigma_clip_level)].index)
#endif

for i in RV_array.index:
    RV_array.loc[i,'Phase'] = ((RV_array.loc[i,'Time'] - Mid_transit)/TESS_period) - int(((RV_array.loc[i,'Time'] - Mid_transit)/TESS_period))
#endfor

unique_days = len(RV_array['Time'].round().unique().tolist())
unique_days_array = RV_array['Time'].round().unique().tolist()

Index3 = np.arange(unique_days)
RV_array_binned_days = pd.DataFrame(index=Index3, columns=['Time', 'Phase', 'RV', 'RV_error'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days):
    RV_array_binned_days.loc[i,'Time'] = np.average(RV_array.Time.drop(RV_array[RV_array.Time.round() != unique_days_array[i]].index))
    
    RV_array_binned_days.loc[i,'Phase'] = ((RV_array_binned_days.loc[i,'Time'] - Mid_transit)/TESS_period) - int(((RV_array_binned_days.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days.loc[i,'RV'] = np.average(RV_array.RV.drop(RV_array[RV_array.Time.round() != unique_days_array[i]].index), 
                                                        weights=1.0/RV_array.RV_error.drop(RV_array[RV_array.Time.round() != unique_days_array[i]].index)**2.0)
                                                        
    RV_array_binned_days.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array.RV_error.drop(RV_array[RV_array.Time.round() != unique_days_array[i]].index)**2.0))/ \
                                                        len(RV_array.RV_error.drop(RV_array[RV_array.Time.round() != unique_days_array[i]].index))
#endfor

#Create plots and output text files
unbinned_RV_plot = directory + planet_name + '_RVplot_unbinned.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array["Time"], RV_array["RV"] - np.median(RV_array.RV), yerr=RV_array["RV_error"], fmt='v',
             ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='SONG')
plt.xlabel('Time (HJD)', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('RVs of ' + planet_name, fontsize=16)
plt.legend()
plt.savefig(unbinned_RV_plot, dpi=150)
plt.close(fig)

unbinned_RV_phase_plot = directory + planet_name + '_RVplot_unbinned_phased.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array["Phase"], RV_array["RV"] - np.median(RV_array.RV), yerr=RV_array["RV_error"], fmt='v',
             ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='SONG')
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=14)
plt.legend()
plt.savefig(unbinned_RV_phase_plot, dpi=150)
plt.close(fig)

Binned_day_RV_plot = directory + planet_name + '_RVplot_binned_day.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_days["Time"], RV_array_binned_days["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days["RV_error"], fmt='v',
             ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='SONG')
plt.xlabel('Time (HJD)', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Day Binned RVs of ' + planet_name, fontsize=16)
plt.legend()
plt.savefig(Binned_day_RV_plot, dpi=150)
plt.close(fig)

Binned_day_RV_phase_plot = directory + planet_name + '_RVplot_binned_day_phased.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_days["Phase"], RV_array_binned_days["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days["RV_error"], fmt='v',
                 ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='Telescope 1')
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Day Binned RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=14)
plt.legend()
plt.savefig(Binned_day_RV_phase_plot, dpi=150)
plt.close(fig)



radvel_out = directory + planet_name + '.rv.csv'
RV_array.to_csv(radvel_out, index=False, columns=['Time', 'RV', 'RV_error', 'Telescope'], header=['time', 'mnvel', 'errvel', 'tel'], na_rep='NaN')

radvel_out = directory + planet_name + '.rv_phase.csv'
RV_array.to_csv(radvel_out, index=False, columns=['Phase', 'RV', 'RV_error', 'Telescope'], header=['phase', 'mnvel', 'errvel', 'tel'], na_rep='NaN')


EXOFAST_RVs_out = directory + planet_name + '.SONG.rv'
EXOFAST_out = open(EXOFAST_RVs_out, 'w')
for i in RV_array.index:
    EXOFAST_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array.loc[i, 'Time'],
                                                                     RV_array.loc[i, 'RV'],
                                                                     RV_array.loc[i, 'RV_error']))
#endfor