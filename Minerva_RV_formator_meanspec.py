import pandas as pd
import pyexcel
import numpy as np
import csv
import matplotlib.pyplot as plt
from numpy import mean

directory = '/Users/u8010412/Desktop/MINERVADATA/MINERVA_Results_July_2019/AU_Mic/'
planet_name = 'AU-Mic-R-M'
input_file = 'AU_Mic_Results_28_06_19.txt'
TESS_period = 8.463213
Mid_transit = 1330.39153
sigma_clip = False
sigma_clip_level = 5.0
Exofast_format = True
Radvel_format = True
Exosam_format = True

Mid_transit = Mid_transit + 2457000.0

raw_data = np.genfromtxt(directory + input_file, dtype='str')
data_length = len(raw_data)

Index = np.arange(data_length)
RV_array = pd.DataFrame(index=Index, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])

RV_array['Time'] = raw_data[:,0]
RV_array['RV'] = raw_data[:,1]
RV_array['RV_error'] = raw_data[:,2]
RV_array['Telescope'] = raw_data[:,3]

RV_array['Time'] = RV_array['Time'].astype(float)
RV_array['RV'] = RV_array['RV'].astype(float)
RV_array['RV_error'] = RV_array['RV_error'].astype(float)
RV_array['Telescope'] = RV_array['Telescope'].astype(float).astype(int)

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

RV_array_telescope1 = RV_array.drop(RV_array[RV_array.Telescope != 1].index)
RV_array_telescope2 = RV_array.drop(RV_array[RV_array.Telescope != 2].index)
RV_array_telescope3 = RV_array.drop(RV_array[RV_array.Telescope != 3].index)
RV_array_telescope4 = RV_array.drop(RV_array[RV_array.Telescope != 4].index)
RV_array_telescope5 = RV_array.drop(RV_array[RV_array.Telescope != 5].index)
RV_array_telescope6 = RV_array.drop(RV_array[RV_array.Telescope != 6].index)

unique_time = len(RV_array['Time'].unique().tolist())
unique_time_array = RV_array['Time'].unique().tolist()
unique_days = len(RV_array['Time'].round().unique().tolist())
unique_days_array = RV_array['Time'].round().unique().tolist()

Index2 = np.arange(unique_time)
RV_array_binned_telescopes = pd.DataFrame(index=Index2, columns=['Time', 'Phase', 'RV', 'RV_error'])
#Bin RV data together from the different telescopes
for i in range(unique_time):
    RV_array_binned_telescopes.loc[i,'Time'] = unique_time_array[i]
    
    RV_array_binned_telescopes.loc[i,'Phase'] = ((unique_time_array[i] - Mid_transit)/TESS_period) - int(((unique_time_array[i] - Mid_transit)/TESS_period))
    
    RV_array_binned_telescopes.loc[i,'RV'] = np.average(RV_array.RV.drop(RV_array[RV_array.Time != unique_time_array[i]].index), 
                                                        weights=1.0/RV_array.RV_error.drop(RV_array[RV_array.Time != unique_time_array[i]].index)**2.0)
                                                        
    RV_array_binned_telescopes.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array.RV_error.drop(RV_array[RV_array.Time != unique_time_array[i]].index)**2.0))/ \
                                                        len(RV_array.RV_error.drop(RV_array[RV_array.Time != unique_time_array[i]].index))
#endfor

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

unique_days_telescope1 = len(RV_array_telescope1['Time'].round().unique().tolist())
unique_days_array_telescope1 = RV_array_telescope1['Time'].round().unique().tolist()
Index4 = np.arange(unique_days_telescope1)
RV_array_binned_days_telescope1 = pd.DataFrame(index=Index4, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope1):
    RV_array_binned_days_telescope1.loc[i,'Time'] = np.average(RV_array_telescope1.Time.drop(RV_array_telescope1[RV_array_telescope1.Time.round() \
                                                               != unique_days_array_telescope1[i]].index))
                                                               
    RV_array_binned_days_telescope1.loc[i,'Phase'] = ((RV_array_binned_days_telescope1.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope1.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope1.loc[i,'RV'] = np.average(RV_array_telescope1.RV.drop(RV_array_telescope1[RV_array_telescope1.Time.round() != unique_days_array_telescope1[i]].index), 
                                                             weights=1.0/RV_array_telescope1.RV_error.drop(RV_array_telescope1[RV_array_telescope1.Time.round() \
                                                             != unique_days_array_telescope1[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope1.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope1.RV_error.drop(RV_array_telescope1[RV_array_telescope1.Time.round() \
                                                                != unique_days_array_telescope1[i]].index)**2.0))/ \
                                                                len(RV_array_telescope1.RV_error.drop(RV_array_telescope1[RV_array_telescope1.Time.round() \
                                                                != unique_days_array_telescope1[i]].index))
    RV_array_binned_days_telescope1.loc[i,'Telescope'] = 1
#endfor

unique_days_telescope2 = len(RV_array_telescope2['Time'].round().unique().tolist())
unique_days_array_telescope2 = RV_array_telescope2['Time'].round().unique().tolist()
Index5 = np.arange(unique_days_telescope2)
RV_array_binned_days_telescope2 = pd.DataFrame(index=Index5, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope2):
    RV_array_binned_days_telescope2.loc[i,'Time'] = np.average(RV_array_telescope2.Time.drop(RV_array_telescope2[RV_array_telescope2.Time.round() \
                                                               != unique_days_array_telescope2[i]].index))
                                                               
    RV_array_binned_days_telescope2.loc[i,'Phase'] = ((RV_array_binned_days_telescope2.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope2.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope2.loc[i,'RV'] = np.average(RV_array_telescope2.RV.drop(RV_array_telescope2[RV_array_telescope2.Time.round() != unique_days_array_telescope2[i]].index), 
                                                             weights=1.0/RV_array_telescope2.RV_error.drop(RV_array_telescope2[RV_array_telescope2.Time.round() \
                                                             != unique_days_array_telescope2[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope2.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope2.RV_error.drop(RV_array_telescope2[RV_array_telescope2.Time.round() \
                                                                != unique_days_array_telescope2[i]].index)**2.0))/ \
                                                                len(RV_array_telescope2.RV_error.drop(RV_array_telescope2[RV_array_telescope2.Time.round() \
                                                                != unique_days_array_telescope2[i]].index))
    RV_array_binned_days_telescope2.loc[i,'Telescope'] = 2
#endfor

unique_days_telescope3 = len(RV_array_telescope3['Time'].round().unique().tolist())
unique_days_array_telescope3 = RV_array_telescope3['Time'].round().unique().tolist()
Index6 = np.arange(unique_days_telescope3)
RV_array_binned_days_telescope3 = pd.DataFrame(index=Index6, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope3):
    RV_array_binned_days_telescope3.loc[i,'Time'] = np.average(RV_array_telescope3.Time.drop(RV_array_telescope3[RV_array_telescope3.Time.round() \
                                                               != unique_days_array_telescope3[i]].index))
                                                               
    RV_array_binned_days_telescope3.loc[i,'Phase'] = ((RV_array_binned_days_telescope3.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope3.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope3.loc[i,'RV'] = np.average(RV_array_telescope3.RV.drop(RV_array_telescope3[RV_array_telescope3.Time.round() != unique_days_array_telescope3[i]].index), 
                                                             weights=1.0/RV_array_telescope3.RV_error.drop(RV_array_telescope3[RV_array_telescope3.Time.round() \
                                                             != unique_days_array_telescope3[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope3.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope3.RV_error.drop(RV_array_telescope3[RV_array_telescope3.Time.round() \
                                                                != unique_days_array_telescope3[i]].index)**2.0))/ \
                                                                len(RV_array_telescope3.RV_error.drop(RV_array_telescope3[RV_array_telescope3.Time.round() \
                                                                != unique_days_array_telescope3[i]].index))
    RV_array_binned_days_telescope3.loc[i,'Telescope'] = 3
#endfor

unique_days_telescope4 = len(RV_array_telescope4['Time'].round().unique().tolist())
unique_days_array_telescope4 = RV_array_telescope4['Time'].round().unique().tolist()
Index7 = np.arange(unique_days_telescope4)
RV_array_binned_days_telescope4 = pd.DataFrame(index=Index7, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope4):
    RV_array_binned_days_telescope4.loc[i,'Time'] = np.average(RV_array_telescope4.Time.drop(RV_array_telescope4[RV_array_telescope4.Time.round() \
                                                               != unique_days_array_telescope4[i]].index))
                                                               
    RV_array_binned_days_telescope4.loc[i,'Phase'] = ((RV_array_binned_days_telescope4.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope4.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope4.loc[i,'RV'] = np.average(RV_array_telescope4.RV.drop(RV_array_telescope4[RV_array_telescope4.Time.round() != unique_days_array_telescope4[i]].index), 
                                                             weights=1.0/RV_array_telescope4.RV_error.drop(RV_array_telescope4[RV_array_telescope4.Time.round() \
                                                             != unique_days_array_telescope4[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope4.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope4.RV_error.drop(RV_array_telescope4[RV_array_telescope4.Time.round() \
                                                                != unique_days_array_telescope4[i]].index)**2.0))/ \
                                                                len(RV_array_telescope4.RV_error.drop(RV_array_telescope4[RV_array_telescope4.Time.round() \
                                                                != unique_days_array_telescope4[i]].index))
    RV_array_binned_days_telescope4.loc[i,'Telescope'] = 4
#endfor

unique_days_telescope5 = len(RV_array_telescope5['Time'].round().unique().tolist())
unique_days_array_telescope5 = RV_array_telescope5['Time'].round().unique().tolist()
Index8 = np.arange(unique_days_telescope5)
RV_array_binned_days_telescope5 = pd.DataFrame(index=Index8, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope5):
    RV_array_binned_days_telescope5.loc[i,'Time'] = np.average(RV_array_telescope5.Time.drop(RV_array_telescope5[RV_array_telescope5.Time.round() \
                                                               != unique_days_array_telescope5[i]].index))
                                                               
    RV_array_binned_days_telescope5.loc[i,'Phase'] = ((RV_array_binned_days_telescope5.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope5.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope5.loc[i,'RV'] = np.average(RV_array_telescope5.RV.drop(RV_array_telescope5[RV_array_telescope5.Time.round() != unique_days_array_telescope5[i]].index), 
                                                             weights=1.0/RV_array_telescope5.RV_error.drop(RV_array_telescope5[RV_array_telescope5.Time.round() \
                                                             != unique_days_array_telescope5[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope5.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope5.RV_error.drop(RV_array_telescope5[RV_array_telescope5.Time.round() \
                                                                != unique_days_array_telescope5[i]].index)**2.0))/ \
                                                                len(RV_array_telescope5.RV_error.drop(RV_array_telescope5[RV_array_telescope5.Time.round() \
                                                                != unique_days_array_telescope5[i]].index))
    RV_array_binned_days_telescope5.loc[i,'Telescope'] = 5
#endfor

unique_days_telescope6 = len(RV_array_telescope6['Time'].round().unique().tolist())
unique_days_array_telescope6 = RV_array_telescope6['Time'].round().unique().tolist()
Index9 = np.arange(unique_days_telescope6)
RV_array_binned_days_telescope6 = pd.DataFrame(index=Index9, columns=['Time', 'Phase', 'RV', 'RV_error', 'Telescope'])
#Create 1 day bins of RV data together from the different telescopes
for i in range(unique_days_telescope6):
    RV_array_binned_days_telescope6.loc[i,'Time'] = np.average(RV_array_telescope6.Time.drop(RV_array_telescope6[RV_array_telescope6.Time.round() \
                                                               != unique_days_array_telescope6[i]].index))
                                                               
    RV_array_binned_days_telescope6.loc[i,'Phase'] = ((RV_array_binned_days_telescope6.loc[i,'Time'] - Mid_transit)/TESS_period) \
                                                      - int(((RV_array_binned_days_telescope6.loc[i,'Time'] - Mid_transit)/TESS_period))
    
    RV_array_binned_days_telescope6.loc[i,'RV'] = np.average(RV_array_telescope6.RV.drop(RV_array_telescope6[RV_array_telescope6.Time.round() != unique_days_array_telescope6[i]].index), 
                                                             weights=1.0/RV_array_telescope6.RV_error.drop(RV_array_telescope6[RV_array_telescope6.Time.round() \
                                                             != unique_days_array_telescope6[i]].index)**2.0)
                                                        
    RV_array_binned_days_telescope6.loc[i,'RV_error'] = np.sqrt(np.sum(RV_array_telescope6.RV_error.drop(RV_array_telescope6[RV_array_telescope6.Time.round() \
                                                                != unique_days_array_telescope6[i]].index)**2.0))/ \
                                                                len(RV_array_telescope6.RV_error.drop(RV_array_telescope6[RV_array_telescope6.Time.round() \
                                                                != unique_days_array_telescope6[i]].index))
    RV_array_binned_days_telescope6.loc[i,'Telescope'] = 6
#endfor

#Create plots and output text files
unbinned_RV_plot = directory + planet_name + '_RVplot_unbinned.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
if len(RV_array_telescope1) > 0:
    plt.errorbar(RV_array_telescope1["Time"], RV_array_telescope1["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope1["RV_error"], fmt='v',
                 ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='Telescope 1')
if len(RV_array_telescope2) > 0:
    plt.errorbar(RV_array_telescope2["Time"], RV_array_telescope2["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope2["RV_error"], fmt='^',
                 ecolor='g', color='g', elinewidth=1, markersize=3, capsize=4, label='Telescope 2')
if len(RV_array_telescope3) > 0:
    plt.errorbar(RV_array_telescope3["Time"], RV_array_telescope3["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope3["RV_error"], fmt='o',
                 ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescope 3')
if len(RV_array_telescope4) > 0:
    plt.errorbar(RV_array_telescope4["Time"], RV_array_telescope4["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope4["RV_error"], fmt='s',
                 ecolor='r', color='r', elinewidth=1, markersize=3, capsize=4, label='Telescope 4')
if len(RV_array_telescope5) > 0:
    plt.errorbar(RV_array_telescope5["Time"], RV_array_telescope5["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope5["RV_error"], fmt='*',
                 ecolor='c', color='c', elinewidth=1, markersize=3, capsize=4, label='Telescope 5')
if len(RV_array_telescope6) > 0:
    plt.errorbar(RV_array_telescope6["Time"], RV_array_telescope6["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope6["RV_error"], fmt='D',
                 ecolor='m', color='m', elinewidth=1, markersize=3, capsize=4, label='Telescope 6')
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
if len(RV_array_telescope1) > 0:
    plt.errorbar(RV_array_telescope1["Phase"], RV_array_telescope1["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope1["RV_error"], fmt='v',
                 ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='Telescope 1')
if len(RV_array_telescope2) > 0:
    plt.errorbar(RV_array_telescope2["Phase"], RV_array_telescope2["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope2["RV_error"], fmt='^',
                 ecolor='g', color='g', elinewidth=1, markersize=3, capsize=4, label='Telescope 2')
if len(RV_array_telescope3) > 0:
    plt.errorbar(RV_array_telescope3["Phase"], RV_array_telescope3["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope3["RV_error"], fmt='o',
                 ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescope 3')
if len(RV_array_telescope4) > 0:
    plt.errorbar(RV_array_telescope4["Phase"], RV_array_telescope4["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope4["RV_error"], fmt='s',
                 ecolor='r', color='r', elinewidth=1, markersize=3, capsize=4, label='Telescope 4')
if len(RV_array_telescope5) > 0:
    plt.errorbar(RV_array_telescope5["Phase"], RV_array_telescope5["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope5["RV_error"], fmt='*',
                 ecolor='c', color='c', elinewidth=1, markersize=3, capsize=4, label='Telescope 5')
if len(RV_array_telescope6) > 0:
    plt.errorbar(RV_array_telescope6["Phase"], RV_array_telescope6["RV"] - np.median(RV_array.RV), yerr=RV_array_telescope6["RV_error"], fmt='D',
                 ecolor='m', color='m', elinewidth=1, markersize=3, capsize=4, label='Telescope 6')
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
if len(RV_array_binned_days_telescope1) > 0:
    plt.errorbar(RV_array_binned_days_telescope1["Time"], RV_array_binned_days_telescope1["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope1["RV_error"], fmt='v',
                 ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='Telescope 1')
if len(RV_array_binned_days_telescope2) > 0:
    plt.errorbar(RV_array_binned_days_telescope2["Time"], RV_array_binned_days_telescope2["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope2["RV_error"], fmt='^',
                 ecolor='g', color='g', elinewidth=1, markersize=3, capsize=4, label='Telescope 2')
if len(RV_array_binned_days_telescope3) > 0:
    plt.errorbar(RV_array_binned_days_telescope3["Time"], RV_array_binned_days_telescope3["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope3["RV_error"], fmt='o',
                 ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescope 3')
if len(RV_array_binned_days_telescope4) > 0:
    plt.errorbar(RV_array_binned_days_telescope4["Time"], RV_array_binned_days_telescope4["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope4["RV_error"], fmt='s',
                 ecolor='r', color='r', elinewidth=1, markersize=3, capsize=4, label='Telescope 4')
if len(RV_array_binned_days_telescope5) > 0:
    plt.errorbar(RV_array_binned_days_telescope5["Time"], RV_array_binned_days_telescope5["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope5["RV_error"], fmt='*',
                 ecolor='c', color='c', elinewidth=1, markersize=3, capsize=4, label='Telescope 5')
if len(RV_array_binned_days_telescope6) > 0:
    plt.errorbar(RV_array_binned_days_telescope6["Time"], RV_array_binned_days_telescope6["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope6["RV_error"], fmt='D',
                 ecolor='m', color='m', elinewidth=1, markersize=3, capsize=4, label='Telescope 6')
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
if len(RV_array_binned_days_telescope1) > 0:
    plt.errorbar(RV_array_binned_days_telescope1["Phase"], RV_array_binned_days_telescope1["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope1["RV_error"], fmt='v',
                 ecolor='b', color='b', elinewidth=1, markersize=3, capsize=4, label='Telescope 1')
if len(RV_array_binned_days_telescope2) > 0:
    plt.errorbar(RV_array_binned_days_telescope2["Phase"], RV_array_binned_days_telescope2["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope2["RV_error"], fmt='^',
                 ecolor='g', color='g', elinewidth=1, markersize=3, capsize=4, label='Telescope 2')
if len(RV_array_binned_days_telescope3) > 0:
    plt.errorbar(RV_array_binned_days_telescope3["Phase"], RV_array_binned_days_telescope3["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope3["RV_error"], fmt='o',
                 ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescope 3')
if len(RV_array_binned_days_telescope4) > 0:
    plt.errorbar(RV_array_binned_days_telescope4["Phase"], RV_array_binned_days_telescope4["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope4["RV_error"], fmt='s',
                 ecolor='r', color='r', elinewidth=1, markersize=3, capsize=4, label='Telescope 4')
if len(RV_array_binned_days_telescope5) > 0:
    plt.errorbar(RV_array_binned_days_telescope5["Phase"], RV_array_binned_days_telescope5["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope5["RV_error"], fmt='*',
                 ecolor='c', color='c', elinewidth=1, markersize=3, capsize=4, label='Telescope 5')
if len(RV_array_binned_days_telescope6) > 0:
    plt.errorbar(RV_array_binned_days_telescope6["Phase"], RV_array_binned_days_telescope6["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days_telescope6["RV_error"], fmt='D',
                 ecolor='m', color='m', elinewidth=1, markersize=3, capsize=4, label='Telescope 6')
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Day Binned RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=14)
plt.legend()
plt.savefig(Binned_day_RV_phase_plot, dpi=150)
plt.close(fig)

Binned_telescopes_RV_plot = directory + planet_name + '_RVplot_binned_telescopes.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_telescopes["Time"], RV_array_binned_telescopes["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_telescopes["RV_error"], fmt='o',
             ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescopes Binned')
plt.xlabel('Time (HJD)', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Binned Telescopes RVs of ' + planet_name, fontsize=16)
plt.savefig(Binned_telescopes_RV_plot, dpi=150)
plt.close(fig)

Binned_telescopes_RV_phase_plot = directory + planet_name + '_RVplot_binned_telescopes_phased.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_telescopes["Phase"], RV_array_binned_telescopes["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_telescopes["RV_error"], fmt='o',
             ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescopes Binned')
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Binned Telescopes RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=14)
plt.savefig(Binned_telescopes_RV_phase_plot, dpi=150)
plt.close(fig)

Binned_telescopes_day_RV_plot = directory + planet_name + '_RVplot_binned_telescopes_day.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_days["Time"], RV_array_binned_days["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days["RV_error"], fmt='o',
             ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescopes and Day Binned')
plt.xlabel('Time (HJD)', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Binned Telescopes \& Day RVs of ' + planet_name, fontsize=16)
plt.savefig(Binned_telescopes_day_RV_plot, dpi=150)
plt.close(fig)

Binned_telescopes_day_RV_phase_plot = directory + planet_name + '_RVplot_binned_telescopes_day_phased.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
plt.errorbar(RV_array_binned_days["Phase"], RV_array_binned_days["RV"] - np.median(RV_array.RV), yerr=RV_array_binned_days["RV_error"], fmt='o',
             ecolor='k', color='k', elinewidth=1, markersize=3, capsize=4, label='Telescopes and Day Binned')
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('Binned Telescopes \& Day RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=12)
plt.savefig(Binned_telescopes_day_RV_phase_plot, dpi=150)
plt.close(fig)

for i in RV_array.index:
    RV_array.loc[i, 'Telescope'] = 'minervaTel' + str(RV_array.loc[i,'Telescope'])
#endfor
#print(RV_array)
radvel_out = directory + planet_name + '.rv.csv'
RV_array.to_csv(radvel_out, index=False, columns=['Time', 'RV', 'RV_error', 'Telescope'], header=['time', 'mnvel', 'errvel', 'tel'], na_rep='NaN')

radvel_out = directory + planet_name + '.rv_phase.csv'
RV_array.to_csv(radvel_out, index=False, columns=['Phase', 'RV', 'RV_error', 'Telescope'], header=['phase', 'mnvel', 'errvel', 'tel'], na_rep='NaN')

radvel_out = directory + planet_name + '.rv_phase_tel_binned.csv'
RV_array_binned_telescopes.to_csv(radvel_out, index=False, columns=['Phase', 'RV', 'RV_error'], header=['phase', 'mnvel', 'errvel'], na_rep='NaN')

radvel_out = directory + planet_name + '.rv_tel_binned.csv'
RV_array_binned_telescopes.to_csv(radvel_out, index=False, columns=['Time', 'RV', 'RV_error'], header=['Time', 'mnvel', 'errvel'], na_rep='NaN')

if len(RV_array_telescope1) > 0:
    EXOFAST_tel1_out_RVs = directory + planet_name + '.MINERVA_Tel1.rv'
    EXOFAST_tel1_out = open(EXOFAST_tel1_out_RVs, 'w')
    for i in RV_array_telescope1.index:
        EXOFAST_tel1_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope1.loc[i, 'Time'],
                                                                                 RV_array_telescope1.loc[i, 'RV'],
                                                                                 RV_array_telescope1.loc[i, 'RV_error']))
    #endfor
#endif

if len(RV_array_telescope2) > 0:
    EXOFAST_tel2_out_RVs = directory + planet_name + '.MINERVA_Tel2.rv'
    EXOFAST_tel2_out = open(EXOFAST_tel2_out_RVs, 'w')
    for i in RV_array_telescope2.index:
        EXOFAST_tel2_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope2.loc[i, 'Time'],
                                                                                 RV_array_telescope2.loc[i, 'RV'],
                                                                                 RV_array_telescope2.loc[i, 'RV_error']))
    #endfor
#endif

if len(RV_array_telescope3) > 0:
    EXOFAST_tel3_out_RVs = directory + planet_name + '.MINERVA_Tel3.rv'
    EXOFAST_tel3_out = open(EXOFAST_tel3_out_RVs, 'w')
    for i in RV_array_telescope3.index:
        EXOFAST_tel3_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope3.loc[i, 'Time'],
                                                                                 RV_array_telescope3.loc[i, 'RV'],
                                                                                 RV_array_telescope3.loc[i, 'RV_error']))
    #endfor
#endif

if len(RV_array_telescope4) > 0:
    EXOFAST_tel4_out_RVs = directory + planet_name + '.MINERVA_Tel4.rv'
    EXOFAST_tel4_out = open(EXOFAST_tel4_out_RVs, 'w')
    for i in RV_array_telescope4.index:
        EXOFAST_tel4_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope4.loc[i, 'Time'],
                                                                                 RV_array_telescope4.loc[i, 'RV'],
                                                                                 RV_array_telescope4.loc[i, 'RV_error']))
    #endfor
#endif

if len(RV_array_telescope5) > 0:
    EXOFAST_tel5_out_RVs = directory + planet_name + '.MINERVA_Tel5.rv'
    EXOFAST_tel5_out = open(EXOFAST_tel5_out_RVs, 'w')
    for i in RV_array_telescope5.index:
        EXOFAST_tel5_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope5.loc[i, 'Time'],
                                                                                 RV_array_telescope5.loc[i, 'RV'],
                                                                                 RV_array_telescope5.loc[i, 'RV_error']))
    #endfor
#endif

if len(RV_array_telescope6) > 0:
    EXOFAST_tel6_out_RVs = directory + planet_name + '.MINERVA_Tel6.rv'
    EXOFAST_tel6_out = open(EXOFAST_tel6_out_RVs, 'w')
    for i in RV_array_telescope6.index:
        EXOFAST_tel6_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array_telescope6.loc[i, 'Time'],
                                                                                 RV_array_telescope6.loc[i, 'RV'],
                                                                                 RV_array_telescope6.loc[i, 'RV_error']))
    #endfor
#endif