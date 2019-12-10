# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:23:38 2018

@author: alexc
"""

import matplotlib.pyplot as plt
import nplab
import scipy.signal as sig
import peakutils
import numpy as np
import datetime
from peakutils.plot import plot as pplot
from parula import cm_data
from matplotlib.colors import LinearSegmentedColormap 
from matplotlib.collections import LineCollection
import os
import Alexsdatafuncs
parula_map = LinearSegmentedColormap.from_list('parula', cm_data)

Alex = Alexsdatafuncs.Alexsdatafuncs()

print """
--------------------------------------------------------
Welcome to the analysis code for HDF5 files produced by
the Ocean Optics Spectrometer
---------------------------------------------------------
"""
raw_input(">> ")

print "\n\n Please select the HDF5 file for analysis"
os.chdir(r'C:\Users\alexc\OneDrive - University Of Cambridge\Documents\PhD Project\Projects')
if nplab.datafile._current_datafile != None:
    nplab.datafile._current_datafile.close()
data_f = nplab.current_datafile(mode = 'r+')
fname = data_f.filename
data_f.show_gui()

print "Would you like to now analyse the data outside of the GUI? y/n"
usr_input = raw_input(">> ")

if usr_input == 'y' or usr_input == 'Y':
    print "OK, would you like to filter noise from the data? y/n"
    usr_input2 = raw_input(">> ")
    print """OK, would you like to output an image with highlighted 
    nanoparticle for each spectrum?"""
    usr_input3 = raw_input(">> ")
else:
    print "Exiting..."
    data_f.close()
    raise(Exception)
    
filelist = data_f['OceanOpticsSpectrometer']
camlist = data_f['LumeneraCamera']

max_peaks = 7
wavelength_start = 420
wavelength_end = 950
peakfind_start = 500
wavelengthvals = filelist[filelist.keys()[0]].attrs['wavelengths']
data_start = (np.abs(wavelengthvals - wavelength_start)).argmin()
data_end = (np.abs(wavelengthvals - wavelength_end)).argmin()
peak_data_start = (np.abs(wavelengthvals - peakfind_start)).argmin()

data_range = data_end - data_start
peak_index_store = [None]*max_peaks
wavelength_store = [None]*data_range
intensity_store = [None]*data_range
peak_wavelengths_store = [None]*max_peaks
peak_intensities_store = [None]*max_peaks

spectrum_number_start = 0
spectrum_number_end = len(filelist)

start_time = datetime.datetime.strptime('2018-01-10T16:26:00.406000', 
                                        "%Y-%m-%dT%H:%M:%S.%f")
finish_time = datetime.datetime.strptime('2018-01-10T17:01:53.517000', 
                                         "%Y-%m-%dT%H:%M:%S.%f")

ref = filelist[filelist.keys()[0]].attrs['reference']
back = filelist[filelist.keys()[0]].attrs['background']
x = wavelengthvals

slash_loc = Alex.slashfinder(fname)
spectrum_indexes = []
fig1=plt.figure(1, figsize = (15,50))
ax1=fig1.add_subplot(111) 

num_peaks_store = []
intensity_initial_peak_store = []
initial_peak_wl_store = []

for j, data_set in enumerate(filelist.values()):
    current_timestamp =  datetime.datetime.strptime(
            data_set.attrs['creation_timestamp'], "%Y-%m-%dT%H:%M:%S.%f" )
    if current_timestamp >= start_time and current_timestamp <= finish_time:
              
        spectrum_indexes = np.append(spectrum_indexes, [j])
        
        #print "j = {0}".format(j)
       
        y_raw = data_set.value
        y = (y_raw-back)/(ref-back)
    
        y_red = y[data_start:data_end]
        x_red = x[data_start:data_end]
        
        #Use Savitsky-Golay noise filter to remove noise if chosen
        
        if usr_input2 == "y" or usr_input2 == "Y":
            y_filtered = sig.savgol_filter(y_red,23,3)
            y_norm = (y_filtered/max(y_filtered))
        else:
            y_norm = (y_red/max(y_red))
            pass
        
        #Set peakfinder to only search above set wavelength value with 
        #peak_data_start
        
    
        peak_indexes = peakutils.indexes(y_norm[peak_data_start:data_end], 
                                         thres = 0.5, min_dist = 30)
        peak_indexes = peak_indexes + peak_data_start
        peak_indexes_orig = peak_indexes
        num_peaks = len(peak_indexes_orig)
    
        if len(peak_indexes) > max_peaks:
            while len(peak_indexes) > max_peaks:
                peak_indexes = np.delete(peak_indexes,[max_peaks])
        elif len(peak_indexes) < max_peaks:
            while len(peak_indexes) < max_peaks:
                peak_indexes = np.append(peak_indexes,[None])
        else:
            pass
        
        peak_wavelengths = []
        for i in range(0,max_peaks):
            if peak_indexes[i] != None:
                peak_wavelengths = np.append(peak_wavelengths, 
                                             [x_red[peak_indexes[i]]])
            else:
                peak_wavelengths = np.append(peak_wavelengths, [None])

        plot_peak_intensities = [None]*max_peaks
        peak_intensities = [None]*max_peaks
        for n in range(0, max_peaks):
            if peak_indexes[n] != None:
                plot_peak_intensities[n] = j + y_norm[peak_indexes[n]]
                peak_intensities[n] = y_norm[peak_indexes[n]]
            else:
                break
        
        #Set gradient linewidths to increase with intensity
        #Set separation of line plots to increase by 1 each time with j
        #
        
        lwidths = y_norm[:-1]
        y_norm_2 = y_norm + j
        points = np.array([x_red, y_norm_2]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis = 1)
        lc = LineCollection(segments, linewidths = lwidths, color = 'blue')
        
        #plot line plots with gradient in line thickness with intensity
        
        ax1.add_collection(lc)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.set_xlim(wavelength_start, wavelength_end)
        ax1.set_ylim(min(spectrum_indexes), max(spectrum_indexes))
        ax1.set_xlabel("Wavelength(nm)", fontsize = 25)
        ax1.set_ylabel("Spectrum number", fontsize = 30)
        ax1.set_title(fname[slash_loc:-3]+"\nNormalised unsorted DF spectra", 
                      fontsize = 30)
            
        #Overlay scatter plot of peak position points
        ax1.scatter(peak_wavelengths, np.full((1,max_peaks), plot_peak_intensities))

        wavelength_store = np.vstack([wavelength_store, x_red])
        intensity_store = np.vstack([intensity_store, y_norm])
        peak_index_store = np.vstack([peak_index_store, peak_indexes])
        peak_wavelengths_store = np.vstack([peak_wavelengths_store, 
                                            peak_wavelengths])
        peak_intensities_store = np.vstack([peak_intensities_store, 
                                            peak_intensities])
        num_peaks_store = np.append(num_peaks_store, num_peaks)
        intensity_initial_peak_store = np.append(intensity_initial_peak_store, 
                                                 peak_intensities[0])
        initial_peak_wl_store = np.append(initial_peak_wl_store,peak_wavelengths[0])
        
        
        if usr_input3 == 'y' or usr_input3 == 'Y':
            picno = 0
            for name in camlist.keys():
                if name == data_set.name[25:]:
                    f, axarr = plt.subplots(2, 1, figsize = (15,15))
                    axarr[1].plot(x_red, y_norm)
                    pplot(x_red, y_norm, peak_indexes_orig)
                    axarr[0].imshow(camlist[camlist.keys()[picno]])
                    axarr[0].set_title(data_set.name[25:], fontsize = 20)
                    axarr[0].text(683, 529, 'O', color='y', fontweight = 'bold')
                    for k in range(0,len(peak_indexes_orig)):
                        axarr[1].text(x_red[peak_indexes_orig[k]], 
                             y_norm[peak_indexes_orig[k]] + 0.005, k)
                        axarr[1].text(max(x_red)+50, max(y_norm)-((k+1)*0.02), 
                                 "Peak " + str(k) + ": " + str(x_red[peak_indexes_orig[k]]) + "nm")     
                    axarr[1].set_ylabel("Normalised intensity")
                    axarr[1].set_xlabel("Wavelength(nm)")    
                else:
                    picno = picno + 1
        else:
            pass
        
    else:
        pass

all_info_wl = [wavelength_store, intensity_store, peak_wavelengths_store, 
               peak_intensities_store]
all_info_indexed = [wavelength_store, intensity_store, peak_index_store, 
                    peak_intensities_store]
sorting_array = np.concatenate((intensity_store, peak_wavelengths_store), axis = 1)
sorted_array = sorted(sorting_array, key=lambda x: x[data_range])

while len(sorted_array[1]) > data_range:
    sorted_array = np.delete(sorted_array, [data_end-data_start], 1)
    
sorted_array = np.delete(sorted_array, 0, 0)

all_sorted_array = np.array(sorted_array, dtype = float)


fig2=plt.figure(2, figsize = (10, 10))
ax2=fig2.add_subplot(111)

mapfig1 = ax2.imshow(all_sorted_array, 
                    extent=[wavelength_start, wavelength_end,
                    0, len(all_sorted_array)], aspect = 'auto', 
                            cmap = parula_map)
    
                        
fig2.colorbar(mapfig1)
ax2.tick_params(axis='x', labelsize=20)
ax2.tick_params(axis='y', labelsize=20)
ax2.set_xlabel("Wavelength(nm)", fontsize = 15)
ax2.set_ylabel("Spectrum number", fontsize = 15)
ax2.set_title(fname[slash_loc:-3] + "\nNormalised DF sorted by peak position", fontsize = 18)

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
    
    fig3 = plt.figure(3, figsize = (10,10))
    ax3 = fig3.add_subplot(111)

    mapfig2 = ax3.imshow(all_sorted_array_2, 
                         extent=[wavelength_start, wavelength_end,
                        0, len(all_sorted_array_2)], aspect = 'auto')
    fig3.colorbar(mapfig2)
    ax3.tick_params(axis='x', labelsize=10)
    ax3.tick_params(axis='y', labelsize=10)
    ax3.set_xlabel("Wavelength(nm)", fontsize = 15)
    ax3.set_ylabel("Spectrum number", fontsize = 15)
    ax3.set_title(fname[slash_loc:-3] + 
                  "\nReduced normalised DF sorted by peak position", 
                  fontsize = 18)
    fig3.savefig(fname[slash_loc:-3] + "_reduced map.png")
else:
    pass

#Calculate more useful information
    
av_num_peaks = np.mean(num_peaks_store)
av_init_peak_pos = np.mean(initial_peak_wl_store)
av_init_peak_intensities = np.mean(intensity_initial_peak_store)

fig1.savefig(fname[slash_loc:-3] + "_normalised data.png")
fig2.savefig(fname[slash_loc:-3] + "_normalised map.png")

data_f.flush()
data_f.close()
