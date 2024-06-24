# -*- coding: utf-8 -*-


import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import pandas as pd
import pyexcel                        # import it to handle CSV files.
import fortranformat as ff
from matplotlib.ticker import AutoMinorLocator

RM_data = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/in_transit_RVs_old.csv'
my_data_plotting_symbols = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/in_transit_symbols.txt'
best_parameters_mcmc_input = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/Temp/best_parameters_mcmc.txt'
residual_array_input = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/residuals.txt'
residual_array_null_input = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/residuals_null.txt'
theoretical_input = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/theoretical_RV.txt'
theoretical_null_input = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/theoretical_RV_null.txt'
plot_out = '/Users/u8010412/Desktop/R-M_stuff/RM_analysis/K2-232_combined_transits/RM_K2-232_combined_transits_plot2.pdf'
lower_xlim = -320
upper_xlim = 350
ymin_value_RV = -45
ymax_value_RV = 40
ymin_value_resid = -35
ymax_value_resid = 50

transit1_symbol = 'o'
transit2_symbol = 'X'
transit3_symbol = '^'

labels_plot = ['29 October 2019', None, '4 January 2020']




def to_bool(value):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))
#enddef

def isint(value):
    try:
        int(value)
        return True
    except:
        return False
#enddef




template_RV_sheet = pyexcel.get_sheet(file_name=RM_data, name_columns_by_row=0)
num_my_rv = pyexcel.Sheet.number_of_rows(template_RV_sheet)
my_datafilelength = num_my_rv
Index = np.arange(num_my_rv)
template_RV = pd.DataFrame(index=Index, columns=['Time', 'RV', 'RV_error'])
template_RV['Time'] = template_RV_sheet.column[0]
template_RV['RV'] = template_RV_sheet.column[1]
template_RV['RV_error'] = template_RV_sheet.column[2]
template_RV = template_RV.astype(float)




my_data_plotting_symbols_array = pd.DataFrame(index=Index, columns=['Symbols'])

if len(my_data_plotting_symbols) >= 1:
    #Read in the plotting symbols.
    f = open(my_data_plotting_symbols, 'r')
    num_my_symbols = 0
    for line in f:
        line = line.strip()
        columns = line.split()
        if isint(columns[0]) == True:
            my_data_plotting_symbols_array.iloc[num_my_symbols] = int(columns[0])
        else:
            my_data_plotting_symbols_array.iloc[num_my_symbols] = columns[0]
        #endelse
        num_my_symbols = num_my_symbols + 1
    #endfor
    f.close()

    if num_my_symbols < my_datafilelength:
        #Less symbols than data points. Fill in the missing symbols using the last symbol in the array.
        for i in range(num_my_symbols, my_datafilelength - 1):
            my_data_plotting_symbols_array.iloc[i] = my_data_plotting_symbols_array.iloc[num_my_symbols - 1]
        # endfor
    #endif
#endif
else:
    my_data_plotting_symbols_array[:] = 16
#endelse


#Read in the best over parameters values.
best_overall_parameters = open(best_parameters_mcmc_input, 'r')

format_read1 = ff.FortranRecordReader('(92X, F15.6)')
line = best_overall_parameters.readline()
min_chi_squared_total = format_read1.read(line)[0]
format_read2 = ff.FortranRecordReader('(92X, I10, 1X, I10)')
format_read3 = ff.FortranRecordReader('(92X, I10)')
line = best_overall_parameters.readline()
loc_min_chi_squared_total = format_read2.read(line)
line = best_overall_parameters.readline()
min_r_chi_squared_total = format_read1.read(line)[0]
line = best_overall_parameters.readline()
loc_min_r_chi_squared_total = format_read2.read(line)
line = best_overall_parameters.readline()
max_likelihood_total = format_read1.read(line)[0]
line = best_overall_parameters.readline()
loc_max_likelihood_total = format_read2.read(line)

line = best_overall_parameters.readline()
best_vsini_min_chi_squared = format_read1.read(line)[0]
line = best_overall_parameters.readline()
best_vsini_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
vsini_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_spin_orbit_min_chi_squared = format_read1.read(line)[0]
line = best_overall_parameters.readline()
best_spin_orbit_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
spin_orbit_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_RV_offset_datasets_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
RV_offset_datasets_stand_dev_mcmc_mean = format_read1.read(line)[0]

line = best_overall_parameters.readline()
best_RV_zero_offset_mcmc_mean = format_read1.read(line)[0]
line = best_overall_parameters.readline()
RV_zero_offset_stand_dev_mcmc_mean = format_read1.read(line)[0]

best_overall_parameters.close()


template_RV = template_RV.sort_values('Time', ascending=True)

my_data_plotting_symbols_array = my_data_plotting_symbols_array.iloc[list(template_RV.index)]

template_RV = template_RV.reset_index(drop=True)
my_data_plotting_symbols_array = my_data_plotting_symbols_array.reset_index(drop=True)

#Apply velocity offsets to RV data.
RV_my_array_column = ["Time", "RV", "RV_error"]
RV_my_array_index = np.arange(my_datafilelength)
RV_my_array = pd.DataFrame(index=RV_my_array_index, columns=RV_my_array_column)

RV_my_array.loc[:, 'Time'] = template_RV.loc[:, 'Time']
RV_my_array.loc[:, 'RV'] = template_RV.loc[:, 'RV'] + best_RV_zero_offset_mcmc_mean + best_RV_offset_datasets_mcmc_mean
RV_my_array.loc[:, 'RV_error'] = template_RV.loc[:, 'RV_error']




columns_for_residual_arrays = ["Time", "RV", "RV_error"]
residual_data = np.genfromtxt(residual_array_input, dtype='float', skip_header=1)
length_residual_data = len(residual_data)
index_residual = np.arange(length_residual_data)
residual_array = pd.DataFrame(index=index_residual, columns=columns_for_residual_arrays, dtype='float')
residual_array['Time'] = residual_data[:,0]
residual_array['RV'] = residual_data[:,1]
residual_array['RV_error'] = residual_data[:,2]


residual_data_null = np.genfromtxt(residual_array_null_input, dtype='float', skip_header=1)
length_residual_data_null = len(residual_data_null)
index_residual_null = np.arange(length_residual_data_null)
residual_array_null = pd.DataFrame(index=index_residual_null, columns=columns_for_residual_arrays, dtype='float')
residual_array_null['Time'] = residual_data_null[:, 0]
residual_array_null['RV'] = residual_data_null[:, 1]
residual_array_null['RV_error'] = residual_data_null[:, 2]


columns_for_RV_array = ["Time", "RV"]
RV_array_data = np.genfromtxt(theoretical_input, dtype='float')
length_RV_array_data = len(RV_array_data)
index_RV_array = np.arange(length_RV_array_data)
RV_array = pd.DataFrame(index=index_RV_array, columns=columns_for_RV_array, dtype='float')
RV_array['Time'] = RV_array_data[:,0]
RV_array['RV'] = RV_array_data[:,1]


RV_null_data = np.genfromtxt(theoretical_null_input, dtype='float')
length_RV_array_null_data = len(RV_null_data)
index_RV_array_null = np.arange(length_RV_array_null_data)
RV_array_null = pd.DataFrame(index=index_RV_array_null, columns=columns_for_RV_array, dtype='float')
RV_array_null['Time'] = RV_null_data[:,0]
RV_array_null['RV'] = RV_null_data[:,1]





RVcolo = ['#ebaf02','#ebaf02',
          '#e45311','#e45311',
          '#ce003d','#ce003d',
          '#9f166a','#9f166a',
          '#644aa0','#644aa0',
          '#148097','#148097',
          '#31b15f','#31b15f',]



unique_symbols = len(my_data_plotting_symbols_array['Symbols'].unique().tolist())
unique_symbols_array = my_data_plotting_symbols_array['Symbols'].unique().tolist()

print(RV_my_array)
print(residual_array)
print(residual_array_null)
print(unique_symbols_array)


#Now just plot the fit with just your RV data.

fig = plt.figure()
gs = gridspec.GridSpec(100, 100)
plt.rc('text', usetex=True)
plt.rc('font', family='Times')
ax1 = plt.subplot(gs[:65, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_RV, ymax_value_RV)
plt.setp(plt.gca(), 'xticklabels', [])
plt.ylabel(r'Radial Velocity (ms$^{-1}$)', fontsize=14)
for i in range(unique_symbols):
    p = 0
    for j in RV_my_array.index:
        if my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit1_symbol:
            color_point = RVcolo[12]
            color_error_point = RVcolo[13]
            labels = labels_plot[0]
            #print(labels)
            #print(p)
        elif my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit2_symbol:
            color_point = RVcolo[1]
            color_error_point = RVcolo[2]
            labels = labels_plot[1]
            #print(labels)
            #print(p)
        elif my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit3_symbol:
            color_point = RVcolo[6]
            color_error_point = RVcolo[7]
            labels = labels_plot[2]
            #print(labels)
            #print(p)
        else:
            color_point = 'k'
            color_error_point = 'k'
            labels = None
            #print(labels)
            #print(p)
    
        if my_data_plotting_symbols_array.loc[j, 'Symbols'] == unique_symbols_array[i]:
            #print(labels)
            #print(p)
            #print(my_data_plotting_symbols_array.loc[j, 'Symbols'])
    
            ax1.errorbar(RV_my_array.loc[j, "Time"]*(24.0*60.0), RV_my_array.loc[j, "RV"], yerr=RV_my_array.loc[j, "RV_error"],
                 fmt=my_data_plotting_symbols_array.loc[j, 'Symbols'], ecolor=color_error_point, color=color_point, mfc=color_point, elinewidth=1.5,
                 markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor=color_point,
                 zorder=1, label=labels if p == 0 else "_nolegend_")
                 
            p = p + 1
#endfor
ax1.plot(RV_array["Time"], RV_array["RV"], linestyle='dashed', linewidth=1.25, color="r", zorder=2)
ax1.plot(RV_array_null["Time"], RV_array_null["RV"], linestyle='-.', linewidth=1.25, color="b",
         zorder=3, alpha=0.6)
ax1.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax1.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax1.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax1.xaxis.set_major_locator(plt.MaxNLocator('auto'))
plt.legend()

ax2 = plt.subplot(gs[65:100, 1:100])
plt.xlim(lower_xlim, upper_xlim)
plt.ylim(ymin_value_resid, ymax_value_resid)
for i in range(unique_symbols):
    p = 0
    for j in RV_my_array.index:
        if my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit1_symbol:
            color_point = RVcolo[12]
            color_error_point = RVcolo[13]
        elif my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit2_symbol:
            color_point = RVcolo[1]
            color_error_point = RVcolo[2]
        elif my_data_plotting_symbols_array.loc[j, 'Symbols'] == transit3_symbol:
            color_point = RVcolo[6]
            color_error_point = RVcolo[7]    
        else:
            color_point = 'k'
            color_error_point = 'k'
            
        if my_data_plotting_symbols_array.loc[j, 'Symbols'] == unique_symbols_array[i]:
            
            ax2.errorbar(residual_array.loc[j, "Time"] * (24.0 * 60.0), residual_array.loc[j, "RV"],
                         yerr=residual_array.loc[j, "RV_error"],
                         fmt=my_data_plotting_symbols_array.loc[j, 'Symbols'], ecolor=color_error_point, color=color_point, mfc=color_point, elinewidth=1.5,
                         markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor=color_point,
                         zorder=1)

            ax2.errorbar(residual_array_null.loc[j, "Time"] * (24.0 * 60.0),
                         residual_array_null.loc[j, "RV"], yerr=residual_array_null.loc[j, "RV_error"],
                         fmt=my_data_plotting_symbols_array.loc[j, 'Symbols'], ecolor='k', color='k', mfc='k',
                         elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.25,
                         markeredgecolor='k', zorder=2)
            p = p + 1

# for i in range(num_my_rv):
#     if my_data_plotting_symbols_array.loc[i, 'Symbols'] == transit1_symbol:
#         color_point = RVcolo[12]
#         color_error_point = RVcolo[13]
#     elif my_data_plotting_symbols_array.loc[i, 'Symbols'] == transit2_symbol:
#         color_point = RVcolo[1]
#         color_error_point = RVcolo[2]
#     elif my_data_plotting_symbols_array.loc[i, 'Symbols'] == transit3_symbol:
#         color_point = RVcolo[6]
#         color_error_point = RVcolo[7]
#     else:
#         color_point = 'k'
#         color_error_point = 'k'

    # ax2.errorbar(residual_array.loc[i, "Time"] * (24.0 * 60.0), residual_array.loc[i, "RV"],
    #              yerr=residual_array.loc[i, "RV_error"],
    #              fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor=color_error_point, color=color_point, mfc=color_point, elinewidth=1.5,
    #              markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.9, markeredgecolor=color_point,
    #              zorder=1)
    #
    # ax2.errorbar(residual_array_null.loc[i, "Time"] * (24.0 * 60.0),
    #              residual_array_null.loc[i, "RV"], yerr=residual_array_null.loc[i, "RV_error"],
    #              fmt=my_data_plotting_symbols_array.loc[i, 'Symbols'], ecolor='k', color='k', mfc='k',
    #              elinewidth=1.5, markersize=5, capsize=0, markeredgewidth=1.0, fillstyle='full', alpha=0.25,
    #              markeredgecolor='k', zorder=2)
    

#endfor
ax2.axhline(y=0, linestyle='dashed', linewidth=1.25, color="r", zorder=2)
plt.xlabel('Time (Minutes From Mid Transit)', fontsize=14)
ax2.tick_params(axis='both', which='both', labelsize=14, direction='in', top="on", right="on", length=6)
ax2.tick_params(which='minor', length=3, color='k', direction='in', top="on", right="on")
ax2.tick_params(which='major', length=6, color='k', direction='in', top="on", right="on", width=1.3)
plt.ylabel(r'O -- C (ms$^{-1}$)', fontsize=14)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.yaxis.set_major_locator(plt.MaxNLocator('auto'))
ax2.xaxis.set_major_locator(plt.MaxNLocator('auto'))
plt.savefig(plot_out, dpi=150, bbox_inches='tight')
plt.close(fig)









# if __name__ == "__main__":
#
#     RVcolo = ['#ebaf02','#ebaf02',
#               '#e45311','#e45311',
#               '#ce003d','#ce003d',
#               '#9f166a','#9f166a',
#               '#644aa0','#644aa0',
#               '#148097','#148097',
#               '#31b15f','#31b15f',]
#
#
#     f, axarr = plt.subplots(2, sharex=True,gridspec_kw={'height_ratios': [3, 1]})
#     axarr[0].plot(rvMod[:,0],rvMod[:,1],'--',c='#777777')
#     axarr[0].errorbar(ferosBin[:,0],ferosBin[:,1],yerr=ferosBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[12],ecolor=RVcolo[13],label='FEROS')
#     axarr[0].errorbar(harpsBin[:,0],harpsBin[:,1],yerr=harpsBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[1],ecolor=RVcolo[2],label='HARPS')
#     axarr[0].errorbar(minAusBin[:,0],minAusBin[:,1],yerr=minAusBin[:,2]*1.8,fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[6],ecolor=RVcolo[7],label='MINERVA-Aus')
#     axarr[0].legend(frameon=False,fontsize=10,ncol=3,loc=3)
#     axarr[0].tick_params(which = 'both', direction = 'in',right=True,top=True)
#
#     axarr[1].plot(rvMod[:,0],rvMod[:,0]*0.0,'--',c='#777777')
#     axarr[1].errorbar(ferosBin[:,0],ferosBin[:,1]-modelFeros,yerr=ferosBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[12],ecolor=RVcolo[13],label='FEROS')
#     axarr[1].errorbar(harpsBin[:,0],harpsBin[:,1]-modelHarps,yerr=harpsBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[1],ecolor=RVcolo[2],label='HARPS')
#     axarr[1].errorbar(minAusBin[:,0],minAusBin[:,1]-modelMinAus,yerr=minAusBin[:,2]*1.8,fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[6],ecolor=RVcolo[7],label='MINERVA-Aus')
#     axarr[1].set_xlabel('BJD$_{\mathregular{TDB}}$',fontsize=16)
#     axarr[0].set_ylabel(r'RV $\mathregular{[ms^{-1}]}$',fontsize=16)
#     axarr[1].set_ylabel(r'O-C',fontsize=16)
#
#     axarr[1].set_xlim(2458450,2458800)
#     plt.tick_params(which = 'both', direction = 'in',right=True,top=True)
#     plt.tight_layout()
#     axarr[1].yaxis.set_minor_locator(MultipleLocator(12.5))
#     f.subplots_adjust(hspace=0)
#     plt.savefig(path+'\\binnedunfolded.pdf')
#     plt.close()
#
#
#     """
#     BINNED FOLDED DATA
#     """
#
#     rvMod[:,0]          = foldAt(rvMod[:,0],period, T0=T_0)
#     minAusBin[:,0]      = foldAt(minAusBin[:,0],period, T0=T_0)
#     harpsBin[:,0]       = foldAt(harpsBin[:,0],period, T0=T_0)
#     ferosBin[:,0]       = foldAt(ferosBin[:,0],period, T0=T_0)
#
#     rvMod = rvMod[rvMod[:,0].argsort()[::-1]]
#
#     f, axarr = plt.subplots(2, sharex=True,gridspec_kw={'height_ratios': [3, 1]})
#     axarr[0].plot(rvMod[:,0],rvMod[:,1],'--',c='#777777')
#     axarr[0].errorbar(ferosBin[:,0],ferosBin[:,1],yerr=ferosBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[12],ecolor=RVcolo[13],label='FEROS')
#     axarr[0].errorbar(harpsBin[:,0],harpsBin[:,1],yerr=harpsBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[1],ecolor=RVcolo[2],label='HARPS')
#     axarr[0].errorbar(minAusBin[:,0],minAusBin[:,1],yerr=minAusBin[:,2]*1.8,fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[6],ecolor=RVcolo[7],label='MINERVA-Aus')
#     axarr[0].legend(frameon=False,fontsize=10,ncol=3,loc=4)
#     axarr[0].tick_params(which = 'both', direction = 'in',right=True,top=True)
#
#     axarr[1].plot(rvMod[:,0],rvMod[:,0]*0.0,'--',c='#777777')
#     axarr[1].errorbar(ferosBin[:,0],ferosBin[:,1]-modelFeros,yerr=ferosBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[12],ecolor=RVcolo[13],label='FEROS')
#     axarr[1].errorbar(harpsBin[:,0],harpsBin[:,1]-modelHarps,yerr=harpsBin[:,2],fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[1],ecolor=RVcolo[2],label='HARPS')
#     axarr[1].errorbar(minAusBin[:,0],minAusBin[:,1]-modelMinAus,yerr=minAusBin[:,2]*1.8,fmt='o',
#                  alpha=alph,markersize=marksize,c=RVcolo[6],ecolor=RVcolo[7],label='MINERVA-Aus')
#     axarr[1].set_xlabel(r'Orbital phase',fontsize=16)
#     axarr[0].set_ylabel(r'RV $\mathregular{[ms^{-1}]}$',fontsize=16)
#     axarr[1].set_ylabel(r'O-C',fontsize=16)
#
#     axarr[1].set_xlim(0,1)
#     plt.tick_params(which = 'both', direction = 'in',right=True,top=True)
#     plt.tight_layout()
#     plt.gca().xaxis.set_minor_locator(MultipleLocator(0.1))
#     axarr[1].yaxis.set_minor_locator(MultipleLocator(12.5))
#     f.subplots_adjust(hspace=0)
#     plt.savefig(path+'\\binnedfolded.pdf')
#     plt.close()
