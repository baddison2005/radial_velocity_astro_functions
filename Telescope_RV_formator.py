import pandas as pd
import pyexcel
import numpy as np
import csv
import matplotlib.pyplot as plt
from numpy import mean

directory = '/Users/u8010412/Desktop/MINERVADATA/MINERVA_Results_July_2019/TOI778/'
planet_name = 'TOI778'
input_file = 'TOI778_TRES.txt'
TESS_period = 4.633751
Mid_transit = 1569.448340
bisector_cut = False
sigma_clip = False
sigma_clip_level = 5.0
change_date_format = False
date_format = 2400000.0
telescope_column = False
bisector_column = False
telescope = 'TRES'
header_in_file = True
#Are the RVs in km/s instead of m/s? Convert them to m/s if in km/s.
RV_format_kms = False
Exofast_format = True
Radvel_format = True
Exosam_format = False
Allesfitter_format = True

Mid_transit = Mid_transit + 2457000.0

if header_in_file:
    raw_data = np.genfromtxt(directory + input_file, dtype='str', comments='#')
    data_length = len(raw_data)
else:
    raw_data = np.genfromtxt(directory + input_file, dtype='str')
    data_length = len(raw_data)

Index = np.arange(data_length)
RV_array = pd.DataFrame(index=Index, columns=['Time', 'Phase', 'RV', 'RV_error', 'Bisector', 'Bisector_error', 'Telescope'])

RV_array['Time'] = raw_data[:,0]
RV_array['RV'] = raw_data[:,1]
RV_array['RV_error'] = raw_data[:,2]
if bisector_column:
    RV_array['Bisector'] = raw_data[:,3]
    RV_array['Bisector_error'] = raw_data[:,4]
else:
    RV_array['Bisector'] = np.nan
    RV_array['Bisector_error'] = np.nan
if telescope_column:
    RV_array['Telescope'] = raw_data[:,5]
else:
    RV_array['Telescope'] = telescope

if change_date_format:
    RV_array['Time'] = RV_array['Time'].astype(float) + date_format
else:
    RV_array['Time'] = RV_array['Time'].astype(float)
if RV_format_kms == True:
    RV_array['RV'] = RV_array['RV'].astype(float)*1000.0
    RV_array['RV_error'] = RV_array['RV_error'].astype(float)*1000.0
else:
    RV_array['RV'] = RV_array['RV'].astype(float)
    RV_array['RV_error'] = RV_array['RV_error'].astype(float)
RV_array['Bisector'] = RV_array['Bisector'].astype(float)
RV_array['Bisector_error'] = RV_array['Bisector_error'].astype(float)
RV_array['Telescope'] = RV_array['Telescope'].astype(str)

#Remove the bad RVs. Bisectors give clue to bad RVs. Use them to remove bad values before sigma clipping.
if bisector_cut == True:
    RV_array = RV_array.drop(RV_array[RV_array.Bisector <= -1.0].index)
    RV_array = RV_array.drop(RV_array[RV_array.Bisector >= 1.0].index)
    RV_array = RV_array.drop(RV_array[RV_array.Bisector.isna()].index)
#endif

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

unique_telescope = len(RV_array['Telescope'].unique().tolist())
unique_telescope_array = RV_array['Telescope'].unique().tolist()
Index_plot = np.arange(10)
plot_array = pd.DataFrame(index=Index_plot, columns=['Symbols', 'Color'])
plot_array['Symbols'] = ['v', '^', 'o', 's', '*', 'D', 'P', 'X', '<', '>']
plot_array['Color'] = ['b', 'g', 'k', 'r', 'c', 'm', 'darkblue', 'darkgreen', 'orange', 'teal']


#Create plots and output text files
unbinned_RV_plot = directory + planet_name + '_RVplot_unbinned.pdf'
fig = plt.figure()
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
for i in range(unique_telescope):
    p = 0
    for j in RV_array.index:
        if RV_array.loc[j, 'Telescope'] == unique_telescope_array[i]:
            plt.errorbar(RV_array.loc[j, "Time"], RV_array.loc[j, "RV"] - np.median(RV_array.RV), yerr=RV_array.loc[j, "RV_error"], fmt=plot_array.loc[i, 'Symbols'],
                         ecolor=plot_array.loc[i, 'Color'], color=plot_array.loc[i, 'Color'], elinewidth=1, markersize=3, capsize=4, label=RV_array.loc[j, "Telescope"] if p == 0 else "_nolegend_")
            p = p + 1
        #endif
    #endfor
#endfor
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
for i in range(unique_telescope):
    p = 0
    for j in RV_array.index:
        if RV_array.loc[j, 'Telescope'] == unique_telescope_array[i]:
            plt.errorbar(RV_array.loc[j, "Phase"], RV_array.loc[j, "RV"] - np.median(RV_array.RV), yerr=RV_array.loc[j, "RV_error"], fmt=plot_array.loc[i, 'Symbols'],
                         ecolor=plot_array.loc[i, 'Color'], color=plot_array.loc[i, 'Color'], elinewidth=1, markersize=3, capsize=4, label=RV_array.loc[j, "Telescope"] if p == 0 else "_nolegend_")
            p = p + 1
        #endif
    #endfor
#endfor
plt.xlabel('Orbital Phase', fontsize=14)
plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
plt.title('RVs of ' + planet_name + ', Phased to TESS Period ' + str(TESS_period) + 'd', fontsize=14)
plt.legend()
plt.savefig(unbinned_RV_phase_plot, dpi=150)
plt.close(fig)

if bisector_column:
    Bisector_plot = directory + planet_name + '_Bisector_span.pdf'
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.rc('font', family='Times')
    for i in range(unique_telescope):
        p = 0
        for j in RV_array.index:
            if RV_array.loc[j, 'Telescope'] == unique_telescope_array[i]:
                plt.errorbar(RV_array.loc[j, "RV"] - np.median(RV_array.RV), RV_array.loc[j, "Bisector"], xerr=RV_array.loc[j, "RV_error"], yerr=RV_array.loc[j, "Bisector_error"], 
                             fmt=plot_array.loc[i, 'Symbols'], ecolor=plot_array.loc[i, 'Color'], color=plot_array.loc[i, 'Color'], elinewidth=1, markersize=3, capsize=4, 
                             label=RV_array.loc[j, "Telescope"] if p == 0 else "_nolegend_")
                p = p + 1
            #endif
        #endfor
    #endfor
    plt.xlabel('Radial Velocity (ms$^{-1}$)', fontsize=14)
    plt.tick_params(axis='both', which='major', direction='in', bottom=True, top=True, left=True, right=True, labelsize=14)
    plt.ylabel(r'Bisector (kms$^{-1}$)', fontsize=14)
    plt.title('Bisector Span of ' + planet_name, fontsize=16)
    plt.legend()
    plt.savefig(Bisector_plot, dpi=150)
    plt.close(fig)
#endif

if Radvel_format:
    radvel_out = directory + planet_name + '.rv.csv'
    RV_array.to_csv(radvel_out, index=False, columns=['Time', 'RV', 'RV_error', 'Telescope'], header=['time', 'mnvel', 'errvel', 'tel'], na_rep='NaN')

    radvel_out = directory + planet_name + '.rv_phase.csv'
    RV_array.to_csv(radvel_out, index=False, columns=['Phase', 'RV', 'RV_error', 'Telescope'], header=['phase', 'mnvel', 'errvel', 'tel'], na_rep='NaN')

if Exofast_format:
    for i in range(unique_telescope):
        EXOFAST_out_RVs = directory + planet_name + '.' + unique_telescope_array[i] + '.rv'
        EXOFAST_out = open(EXOFAST_out_RVs, 'w')
        for j in RV_array.index:
            if RV_array.loc[j, 'Telescope'] == unique_telescope_array[i]:
                EXOFAST_out.write('{:<13.6f}     {:<8.1f}     {:<5.1f}\n'.format(RV_array.loc[j, 'Time'],
                                                                                 RV_array.loc[j, 'RV'],
                                                                                 RV_array.loc[j, 'RV_error']))
            #endif
        #endfor
    #endfor
#endif

if Allesfitter_format:
    for i in range(unique_telescope):
        allesfitter_out = directory + unique_telescope_array[i] + '.csv'
        RV_array['RV'] = RV_array['RV']/1000.0
        RV_array['RV_error'] = RV_array['RV_error']/1000.0
        RV_array.to_csv(allesfitter_out, index=False, columns=['Time', 'RV', 'RV_error'], header=['# time', 'rv', 'rv_err'], na_rep='NaN')