# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:23:38 2018

@author: alexc
"""


import matplotlib.pyplot as plt
import nplab
import scipy.signal as sig
import scipy.constants as const
import peakutils
import numpy as np
import datetime
from peakutils.plot import plot as pplot
from parula import cm_data
from matplotlib.colors import LinearSegmentedColormap 
from matplotlib.collections import LineCollection
import os
import Alexsdatafuncs
import matplotlib.cm as cm
import time
parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
copper_map = LinearSegmentedColormap.from_list('copper', cm_data)

Alex_data = Alexsdatafuncs.data_manipulators()
Alex_plot = Alexsdatafuncs.plotting_funcs()



print """
--------------------------------------------------------
Welcome to the general analysis code for HDF5 files produced by
the Ocean Optics Spectrometer
---------------------------------------------------------
"""
raw_input("Press any key to select HDF5 file >> ")

os.chdir(r"""C:\Users\alexc\OneDrive - University Of Cambridge\Documents\PhD Project\Projects""")

if nplab.datafile._current_datafile != None:
    nplab.datafile._current_datafile.close()
data_f = nplab.current_datafile(mode = 'r')
fname = data_f.filename

questions = ["Would you like to use default plot settings?",
        """Should we plot spectra taken over a range of time ('y') or with
             certain names ('n')?""",
             "OK, would you like to filter noise from the data? y/n",
             """OK, would you like to output an image with highlighted 
             nanoparticle for each spectrum?""",
             "Would you like to complete a gaussian fit on your data?",
             "Would you like to plot you data with line gradients?",
             "Would you like to plot the peak positions too?"]

answers = [True]*len(questions)

def ans_to_bool(ans):
    if ans in ["y","Y"]:
        return True
    elif ans in ["n", "N"]:
        return False
    elif ans in ["e", "E"]:
        raise ValueError("Exiting...")
        data_f.close()
        quit()     

for i,q in enumerate(questions):
    print q
    answers[i]=ans_to_bool(raw_input(">> "))
    if answers[0] == True:
        print "Using default settings...."
        answers = [False, True, True, False, False, False, False]
        break
        

filelist = data_f['OceanOpticsSpectrometer']
#camlist = data_f['LumeneraCamera']

max_peaks = 7
wavelength_start = 405
wavelength_end = 1000
peakfind_start = 500
wavelengthvals = filelist[filelist.keys()[0]].attrs['wavelengths']
data_start = (np.abs(wavelengthvals - wavelength_start)).argmin()
data_end = (np.abs(wavelengthvals - wavelength_end)).argmin()
peak_data_start = (np.abs(wavelengthvals - peakfind_start)).argmin()

data_range = data_end - data_start
peak_index_store = []
wavelength_store = []
name_store = []
intensity_store_norm = []
intensity_store = []
peak_wavelengths_store = []
peak_intensities_store = []


if answers[1] == True:
    start_time = datetime.datetime.strptime('2017-12-05T15:51:36.805000', 
                                        "%Y-%m-%dT%H:%M:%S.%f")
    finish_time = datetime.datetime.strptime('2017-12-05T16:02:02.837000', 
                                        "%Y-%m-%dT%H:%M:%S.%f") 
    search_str = "_"    
elif answers[1] == False:
    time_list = [] 
    name_list =[]
    for i in range(0,len(filelist)):
        time_i = datetime.datetime.strptime(
                filelist.values()[i].attrs['creation_timestamp'][:19], 
                "%Y-%m-%dT%H:%M:%S")
        time_s = time.mktime(time_i.timetuple())
        time_list = np.append(time_list, time_s)
        name_list = np.append(name_list, time_i)
    start_time = name_list[np.argmin(time_list)]
    finish_time = name_list[np.argmax(time_list)]

    search_str = raw_input("Plot files containing string in title >> ")
    
try:
    for i in range(0, len(filelist)):
        try:
            ref = filelist[filelist.keys()[i]].attrs['reference']
            print "Reference spectrum taken from spectrum: " 
            print filelist.keys()[i]
            break
        except KeyError:
            pass
except KeyError:
    print """Coudn't find a stored reference attribute, \n
    find the file containing a reference and paste this into the code"""
    ref = [0]*data_range
    
    
back = filelist[filelist.keys()[0]].attrs['background']
x = wavelengthvals

slash_loc = Alex_data.slashfinder(fname)

spectrum_indexes = []
num_peaks_store = []
intensity_initial_peak_store = []
initial_peak_wl_store = []
second_peak_wl_store = []

peak_width_guess = 50
exp_number_of_peaks = 2

no_spectra=-1

for data_set in filelist.values():
    
    current_timestamp =  datetime.datetime.strptime(
            data_set.attrs['creation_timestamp'][:19], "%Y-%m-%dT%H:%M:%S")
    
    no_spectra += 1
    
    if current_timestamp >= start_time and current_timestamp <= finish_time and (search_str in filelist.keys()[no_spectra]):
        
        print "Current spectrum: " + data_set.name[25:]
        
        spectrum_indexes = np.append(spectrum_indexes, [no_spectra])
       
        y_raw = np.array(data_set.value)
        y = (y_raw-back)/(ref-back)
    
        y_red = y[data_start:data_end]
        x_red = x[data_start:data_end]
        
        #Use Savitsky-Golay noise filter to remove noise if chosen
        
        if answers[2] == True:
            y_filtered = sig.savgol_filter(y_red,23,3)
            y_filtered_corr = y_filtered - min(y_filtered)
            y_norm = (y_filtered_corr/max(y_filtered_corr))
            y_red = y_filtered
        else:
            y_red = y[data_start:data_end] - min(y[data_start:data_end])
            y_norm = (y_red/max(y_red))
            pass
        
        #Set peakfinder to only search above set wavelength value with 
        #peak_data_start
        
        #Automatically search for the correct threshold leading the code to 
        #find the defined number of peaks
        thresh_guess = 0.2
        find_peaks = False
        while find_peaks == False and exp_number_of_peaks != 0:
            peak_indexes = peakutils.indexes(y_norm[peak_data_start:data_end], 
                                             thres = thresh_guess, 
                                             min_dist = peak_width_guess*2)
            if len(peak_indexes) >= exp_number_of_peaks:
                find_peaks = True
            elif thresh_guess <=0.01:
                while find_peaks == False:
                    peak_indexes = peakutils.indexes(y_norm[peak_data_start:data_end], 
                                                     thres = thresh_guess, 
                                                     min_dist = peak_width_guess*2)
                    if len(peak_indexes) >= exp_number_of_peaks:
                        find_peaks = True
                    elif thresh_guess > 1:
                        print "Didn't find threshold for spectrum" + str(no_spectra)
                        print "Need to change expected number or peaks"
                        thresh_guess = 0.5
                        break
                    else:
                        thresh_guess += 0.01
                break
            else:
                thresh_guess -= 0.01
        
        peak_indexes = peak_indexes + peak_data_start
        peak_indexes_orig = peak_indexes
        num_peaks = len(peak_indexes_orig)
    
        if len(peak_indexes) > max_peaks:
            while len(peak_indexes) > max_peaks:
                peak_indexes = np.delete(peak_indexes,[max_peaks])
        elif len(peak_indexes) < max_peaks:
            while len(peak_indexes) < max_peaks:
                peak_indexes = np.append(peak_indexes,[0])
        else:
            pass
        
        peak_wavelengths = []
        for i in range(0,max_peaks):
            if peak_indexes[i] != 0:
                peak_wavelengths = np.append(peak_wavelengths, 
                                             [x_red[peak_indexes[i]]])
            else:
                peak_wavelengths = np.append(peak_wavelengths, [0])
                
        peak_intensities = []        
        for i in range(0, max_peaks):
            if peak_indexes[i] != 0:
                peak_intensities = np.append(peak_intensities, 
                                             [y_norm[peak_indexes[i]]])
            else:
                peak_intensities = np.append(peak_intensities, [0])
         
        #Gaussian fitting    
        if answers[4]==True:
            try:
                multi_gauss_fit_a = Alex_data.gaussian_fit_2(y_norm, x_red, 
                                                             peak_indexes_orig, 
                                                             peak_intensities)
            except UnboundLocalError:
                print "Gauss fit didnt work for spectrum" + str(no_spectra) 
                pass
        else:
            pass
        
        wavelength_store.append(x_red)
        name_store.append(data_set.name[25:])
        intensity_store.append(y_red)
        intensity_store_norm.append(y_norm)
        peak_index_store.append(peak_indexes)
        peak_wavelengths_store.append(peak_wavelengths)
        peak_intensities_store.append(peak_intensities)
        num_peaks_store = np.append(num_peaks_store, num_peaks)
        intensity_initial_peak_store = np.append(intensity_initial_peak_store, 
                                                 peak_intensities[0])
        initial_peak_wl_store = np.append(initial_peak_wl_store,peak_wavelengths[0])
        second_peak_wl_store = np.append(second_peak_wl_store, peak_wavelengths[1])
        
        
        #Plot each individual spectrum with a picture showing where the 
        #data was taken from
#        picno = 0
#        if answers[3] == True:
#            for name in camlist.keys():
#                if name == data_set.name[25:]:
#                    fig, axarr = plt.subplots(2, 1, figsize = (7,7))
#                    axarr[1].plot(x_red, y_red)
#                    pplot(x_red, y_red, peak_indexes_orig)
#                    axarr[0].imshow(camlist[camlist.keys()[picno]])
#                    axarr[0].set_title(data_set.name[25:], fontsize = 20)
#                    axarr[0].text(683, 529, 'O', color='y', fontweight = 'bold')
#                    for k in range(0,len(peak_indexes_orig)):
#                        axarr[1].text(x_red[peak_indexes_orig[k]], 
#                             y_norm[peak_indexes_orig[k]] + 0.005, k)
#                        axarr[1].text(max(x_red)+50, max(y_norm)-((k+1)*0.1), 
#                                 "Peak " + str(k) + ": " + str(x_red[peak_indexes_orig[k]]) + "nm")     
#                    axarr[1].set_ylabel("Intensity referenced to Au")
#                    axarr[1].set_xlabel("Wavelength(nm)")  
#                    print x_red[peak_indexes_orig[0]]
#                else:
#                    picno = picno + 1
        
#        else:
#            pass
        
    else:
        pass

#Not sure whether this is necessary, but for formatting
        
wavelength_store = np.array(wavelength_store)        
intensity_store = np.array(intensity_store)
name_store = np.array(name_store)
intensity_store_norm = np.array(intensity_store_norm)
peak_index_store = np.array(peak_index_store)
peak_wavelengths_store = np.array(peak_wavelengths_store)
peak_intensities_store = np.array(peak_intensities_store)
num_peaks_store = np.array(num_peaks_store)
intensity_initial_peak_store = np.array(intensity_initial_peak_store)
initial_peak_wl_store = np.array(initial_peak_wl_store)
second_peak_wl_store = np.array(second_peak_wl_store)


#Calculate useful values

av_num_peaks = np.mean(num_peaks_store)
av_second_peak_pos = np.mean(second_peak_wl_store)
av_init_peak_pos_stddev = np.std(initial_peak_wl_store)
av_second_peak_pos_stdev = np.std(second_peak_wl_store)
av_init_peak_pos = 660

init_peak = []
second_peak = []
for i in range(0, len(initial_peak_wl_store)): 
    if initial_peak_wl_store[i] < (av_init_peak_pos + av_init_peak_pos_stddev):
        init_peak = np.append(init_peak, initial_peak_wl_store[i])
    else:
        pass
av_init_peak = np.mean(init_peak)
for i in range(0, len(second_peak_wl_store)): 
    if initial_peak_wl_store[i] < (av_second_peak_pos + av_second_peak_pos_stdev):
        second_peak = np.append(second_peak, second_peak_wl_store[i])
    else:
        pass
    
#Draw line through average    
if answers[6] == True:    
    ax1.plot([av_init_peak]*len(np.linspace(64,70)), np.linspace(64,70),'r--')
    ax1.plot([av_second_peak_pos]*len(np.linspace(64,70)), np.linspace(64,70),'r--')
else:
    pass

av_init_peak_intensities = np.mean(intensity_initial_peak_store)
average_all_intensities_norm = np.mean(intensity_store_norm, axis=0)
average_all_intensities = np.mean(intensity_store, axis=0)
av_stddev_norm = (np.std((intensity_store_norm), axis=0, dtype = np.float64))/np.sqrt(no_spectra)
av_stddev = np.std((intensity_store), axis=0, dtype = np.float64)/np.sqrt(no_spectra)


min_peak_pos_av = x_red[np.argmin(average_all_intensities_norm)]

#Plot the normal of the average of all spectra
picno = 0
fig1=plt.figure(picno+1)
ax1=fig1.add_subplot(111)
tax1=ax1.twiny()
ax1.plot(x_red, average_all_intensities_norm*100, label = 'Averaged all intensities')
ax1.plot(x_red, av_stddev_norm*100 + average_all_intensities_norm*100, 'r--', 
         label = "Stddev")
ax1.plot(x_red, average_all_intensities_norm*100 - av_stddev_norm*100, 'r--')
ax1.fill_between(x_red, average_all_intensities_norm*100 - av_stddev_norm*100,
                av_stddev_norm*100 + average_all_intensities_norm*100, color = "grey", 
                alpha = 0.5, zorder = 3)
ax1.legend(loc = "best")
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set_xlabel("Wavelength(nm)", fontsize = 16)
ax1.set_ylabel("Absolute normalised scattering (%)", fontsize = 16)
ax1.set_title(fname[slash_loc:-3] + "\n Averaged normalised spectra", 
              fontsize = 18, y=1.2)
ax1.set_xlim([min(x_red), max(x_red)])
ax1.set_ylim(0, max(av_stddev_norm*100 + average_all_intensities_norm*100)*1.05)

# Plot the average values of all spectra

fig2=plt.figure(picno+2)
ax2=fig2.add_subplot(111)
tax2 = ax2.twiny()
ax2.plot(x_red, average_all_intensities*100, "green", label = "Average NP scattering")
ax2.plot(x_red, av_stddev*100 + average_all_intensities*100, 'r--', label = "Stddev NP scattering")
ax2.plot(x_red, average_all_intensities*100 - av_stddev*100, 'r--')
ax2.fill_between(x_red, average_all_intensities*100 - av_stddev*100,
                av_stddev*100 + average_all_intensities*100, color = "orange", 
                alpha = 0.3, zorder = 3)
ax2.legend(loc = "best")
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.set_xlabel("Wavelength(nm)", fontsize = 16)
ax2.set_ylabel("Absolute scattering (%)", fontsize = 16)
ax2.set_ylim([0, max(av_stddev*100 + average_all_intensities*100)*1.05])
ax2.set_xlim([min(x_red), max(x_red)])
ax2.set_title(fname[slash_loc:-3] + "\n Averaged spectra", 
              fontsize = 18, y=1.2)
    
#Choose energies to be plotted on top axis

energies = np.array([3, 2.5, 2.25, 2, 1.75, 1.5, 1.25, 1])

#Create array of energies for all wavelengths
Alex_plot.energy_ticks(x_red)

j = 0
new_tick_locations = []    
for i in range(0,len(x_red)):
    if x_red[i] > Alex_plot.e_wl_convert(energies)[j]:
        new_tick_locations = np.append(new_tick_locations, x_red[i])
        j += 1
    else:
        pass

tax2.set_xlim(ax2.get_xlim())
tax2.set_xticks(new_tick_locations)
tax2.set_xticklabels(Alex_plot.energy_ticks(new_tick_locations))
tax2.set_xlabel(r"Energy(eV)")

tax1.set_xlim(ax1.get_xlim())
tax1.set_xticks(new_tick_locations)
tax1.set_xticklabels(Alex_plot.energy_ticks(new_tick_locations))
tax1.set_xlabel(r"Energy(eV)")

#store all information in matrices for access later

all_info_wl = [name_store, wavelength_store, intensity_store, peak_wavelengths_store, 
               peak_intensities_store]
all_info_wl_norm = [name_store, wavelength_store, intensity_store_norm, peak_wavelengths_store, 
               peak_intensities_store]
all_info_indexed_norm = [name_store, wavelength_store, intensity_store_norm, peak_index_store, 
                    peak_intensities_store]

#need to ignore transverse mode for peak sorting, so sort by coupled modes instead

sorted_array = Alex_data.sort_coupled_modes(intensity_store_norm, 
                                            peak_wavelengths_store)
sorted_pk_intensities = Alex_data.sort_coupled_modes(peak_intensities_store, 
                                                     peak_wavelengths_store)

#ax1.scatter(peak_pos, y_value[peak_pos])
sorted_peak_wls = []
for i in range(0, len(sorted_array)):
    sorted_peak_wls.append(sorted_array[i][len(intensity_store_norm):])
    
sorted_peak_wls = np.array(sorted_peak_wls)

while len(sorted_array[1]) > data_range:
    sorted_array = np.delete(sorted_array, [data_range], 1)
    
while len(sorted_pk_intensities[1]) > max_peaks:    
    sorted_pk_intensities = np.delete(sorted_pk_intensities, [max_peaks], 1)    

#Use function to add offset between each spectrum
off_fact = 0.4
sorted_array_offset = Alex_plot.offset_spectra(sorted_array, off_fact)
#sorted_pks_offset = Alex_plot.offset_spectra(sorted_pk_intensities, 0.1)

#Plot sorted and offset spectra with line thickness as a function of intensity

fig3 = plt.figure(picno+3, figsize = (7,0.5*len(sorted_array)))
ax3 = fig3.add_subplot(111)
tax3 = ax3.twiny()
tax3.set_xlim(ax3.get_xlim())
tax3.set_xticks(new_tick_locations)
tax3.set_xticklabels(Alex_plot.energy_ticks(new_tick_locations))
tax3.set_xlabel(r"Energy(eV)")
e_red = 1243.125/x_red
#e_red = np.flip(e_red, axis = 0)


colours = cm.viridis(np.linspace(0,1,len(sorted_array_offset[0])))
#colours = parula_map(np.linspace(0,1,len(sorted_array_offset)))
for i in range(0, len(sorted_array_offset)):
    lwidths = (sorted_array_offset[i] - i*off_fact)*3
    points = np.array([e_red, sorted_array_offset[i]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis = 1)
    lc_2 = LineCollection(segments, linewidths = lwidths, color = colours[600])
    ax3.add_collection(lc_2)
    spec1 = np.array(sorted_array_offset[i], dtype = float)
    spec2 = np.ones(len(sorted_array_offset[i]), dtype = float)*i*off_fact
    ax3.fill_between(e_red, spec1, spec2, color = colours[600], 
                alpha = 0.2, zorder = 3)
    #ax3.scatter(sorted_peak_wls[i], np.full((1,max_peaks), sorted_pks_offset[i]))


   
ax3.set_ylim(0, len(sorted_array_offset)*off_fact + 1 )
ax3.set_xlim(min(e_red), max(e_red))
ax3.set_xlabel("Wavelength (nm)", fontsize = 16)
ax3.set_ylabel("Spectrum_no + Scattered intensity", fontsize = 16)

tax3.set_xlim(ax3.get_xlim())
tax3.set_xticks(new_tick_locations)
tax3.set_xticklabels(Alex_plot.energy_ticks(new_tick_locations))
tax3.set_xlabel(r"Energy(eV)")

#plot map of sorted spectra

fig4=plt.figure(picno+4, figsize = (10, 15))
ax4=fig4.add_subplot(111)

all_sorted_array = np.array(sorted_array, dtype=float)

mapfig1 = ax4.imshow(all_sorted_array, 
                    extent=[wavelength_start, wavelength_end,
                    0, len(all_sorted_array)], aspect = 'auto', 
                            cmap = parula_map, vmin = min(average_all_intensities))
    
                        
fig4.colorbar(mapfig1)
ax4.tick_params(axis='x', labelsize=16)
ax4.tick_params(axis='y', labelsize=16)
ax4.set_xlabel("Wavelength(nm)", fontsize = 16)
ax4.set_ylabel("Spectrum number", fontsize = 16)
ax4.set_title(fname[slash_loc:-3] + "\nNormalised DF sorted by peak position", 
              fontsize = 18, y=1.2)
tax4 = ax4.twiny()
tax4.set_xlim(ax3.get_xlim())
tax4.set_xticks(new_tick_locations)
tax4.set_xticklabels(Alex_plot.energy_ticks(new_tick_locations))
tax4.set_xlabel(r"Energy(eV)")


print """Would you like to remove data above/below a max/min 
    peak value?"""
usr_input4 = raw_input(">> ")
if usr_input4 == 'y' or usr_input4 == 'Y':
    print "Remove data with initial peak wl (nm) below: "
    min_pk_wl = input(">> ")
    print "Remove data with initial peak wl (nm) above: "
    max_pk_wl = input(">> ")

    wavelength_store_cut = wavelength_store
    intensity_store_cut = intensity_store
    peak_wl_store_cut = peak_wavelengths_store

    del_no = 0
    while del_no < len(peak_wl_store_cut):
        if peak_wl_store_cut[del_no][0] < min_pk_wl or peak_wl_store_cut[del_no][0] > max_pk_wl:
            wavelength_store_cut = np.delete(wavelength_store_cut, del_no, 0)
            intensity_store_cut = np.delete(intensity_store_cut, del_no, 0)
            peak_wl_store_cut = np.delete(peak_wl_store_cut, del_no, 0)
        else:
            del_no = del_no + 1
        
    all_info_wl_cut = [wavelength_store_cut, intensity_store_cut, peak_wl_store_cut]

    sorting_array_2 = np.concatenate((intensity_store_cut, peak_wl_store_cut), axis = 1)
    sorted_array_2 = sorted(sorting_array_2, key=lambda x: x[data_range])
    
    while len(sorted_array_2[1]) > data_range:
        sorted_array_2 = np.delete(sorted_array_2, [data_end-data_start], 1)
    
    sorted_array_2 = np.delete(sorted_array_2, 0, 0)
    all_sorted_array_2 = np.array(sorted_array_2, dtype = float)
    
    fig5 = plt.figure(picno+5, figsize = (10,10))
    ax5 = fig5.add_subplot(111)

    mapfig2 = ax5.imshow(all_sorted_array_2, 
                         extent=[wavelength_start, wavelength_end,
                        0, len(all_sorted_array_2)], aspect = 'auto')
    fig5.colorbar(mapfig2)
    ax5.tick_params(axis='x', labelsize=10)
    ax5.tick_params(axis='y', labelsize=10)
    ax5.set_xlabel("Wavelength(nm)", fontsize = 15)
    ax5.set_ylabel("Spectrum number", fontsize = 15)
    ax5.set_title(fname[slash_loc:-3] + 
                 "\nReduced normalised DF sorted by peak position", 
                     fontsize = 18)
    fig5.savefig(fname[slash_loc:-3] + "_reduced map.png")
else:
    pass

data_f.flush()
data_f.close()
