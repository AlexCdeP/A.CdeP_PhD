
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib.pyplot as plt
import nplab
import scipy.signal as sig
import peakutils
import numpy as np
import datetime
from peakutils.plot import plot as pplot
from parula import cm_data
from matplotlib.colors import LinearSegmentedColormap 
import os
parula_map = LinearSegmentedColormap.from_list('parula', cm_data)


print """
--------------------------------------------------------
Welcome to the analysis code for HDF5 files produced by
the Ocean Optics Spectrometer
---------------------------------------------------------
"""

print "\n\n Please select the HDF5 file for analysis"
os.chdir(r'C:\Users\alexc\OneDrive - University Of Cambridge\Documents\PhD Project\Projects')
if nplab.datafile._current_datafile != None:
    nplab.datafile._current_datafile.close()
data_f = nplab.current_datafile(mode = 'r+')
#data_f = nplab.datafile.DataFile(fname,'r+')
fname = data_f.filename
data_f.show_gui()

print "Would you like to now analyse the data outside of the GUI? y/n"

usr_input = raw_input(">> ")

if usr_input == 'y' or usr_input == 'Y':
    print "OK, would you like to filter noise from the data? y/n"
    usr_input2 = raw_input(">> ")
    print """OK, would you like to output an image with highlighted nanoparticle
            for each spectrum?"""
    usr_input3 = raw_input(">> ")
else:
    print "Exiting..."
    data_f.close()
    raise(Exception)
    
filelist = data_f['OceanOpticsSpectrometer']
camlist = data_f['LumeneraCamera']

fig1=plt.figure(1)
ax1=fig1.add_subplot(111)
fig2=plt.figure(2)
ax2=fig2.add_subplot(111)
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)

max_peaks = 7
wavelength_start = 420
wavelength_end = 950
wavelengthvals = filelist[filelist.keys()[0]].attrs['wavelengths']
data_start = (np.abs(wavelengthvals - wavelength_start)).argmin()
data_end = (np.abs(wavelengthvals - wavelength_end)).argmin()

data_range = data_end - data_start
peak_index_store = [None]*max_peaks
wavelength_store = [None]*data_range
intensity_store = [None]*data_range
peak_wavelengths_store = [None]*max_peaks

spectrum_number_start = 1
spectrum_number_end = len(filelist)

start_time = datetime.datetime.strptime('2017-12-13T11:28:48.170000', 
                                        "%Y-%m-%dT%H:%M:%S.%f")
finish_time = datetime.datetime.strptime('2017-12-13T11:47:29.227000', 
                                         "%Y-%m-%dT%H:%M:%S.%f")

ref = filelist[filelist.keys()[0]].attrs['reference']
back = filelist[filelist.keys()[0]].attrs['background']
x = wavelengthvals

for data_set in filelist.values():
    current_timestamp =  datetime.datetime.strptime(
            data_set.attrs['creation_timestamp'], "%Y-%m-%dT%H:%M:%S.%f" )
    if current_timestamp >= start_time and current_timestamp <= finish_time:
       
        picno = 0
        for name in camlist.keys():
            if name == data_set.name[25:]:
                plt.figure(figsize=(10,6))
                plt.title(data_set.name[25:], fontsize = 20)
                plt.text(683, 529, 'O', color='y', fontweight = 'bold')
                plt.imshow(camlist[camlist.keys()[picno]])
            else:
                picno = picno + 1
        
        y_raw = data_set.value 
        y = (y_raw-back)/(ref-back)
    
        y_red = y[data_start:data_end]
        x_red = x[data_start:data_end]

        if usr_input2 == "y" or usr_input2 == "Y":
            y_filtered = sig.savgol_filter(y_red,23,3)
        else:
            pass
    
        y_norm = (y_filtered/max(y_filtered))
    
        peak_indexes = peakutils.indexes(y_norm, thres = 0.5, min_dist = 30)
        peak_indexes_orig = peak_indexes
    
        if len(peak_indexes) > max_peaks:
            while len(peak_indexes) > max_peaks:
                peak_indexes = np.delete(peak_indexes,[max_peaks])
        elif len(peak_indexes) < max_peaks:
            while len(peak_indexes) < max_peaks:
                peak_indexes = np.append(peak_indexes,[None])
        else:
            pass
        
        i=0
        peak_wavelengths = []
        while i < max_peaks:
            peak_wavelengths = np.append(peak_wavelengths, [x_red[peak_indexes[i]]])
            i = i + 1
        
        ax1.plot(x_red, y_norm, 'k', lw=2)
       
        ax1.tick_params(axis='x', labelsize=10)
        ax1.tick_params(axis='y', labelsize=10)
       
        ax1.set_xlim(wavelength_start, wavelength_end)
        ax1.set_ylim(0,1)
        
        ax1.set_xlabel("Wavelength(nm)", fontsize = 10)
        ax1.set_ylabel("Normalised Scattered Intensity", fontsize = 10)
        ax1.set_title("Normalised DF spectra", fontsize = 18)
    
        wavelength_store = np.vstack([wavelength_store, x_red])
        intensity_store = np.vstack([intensity_store, y_norm])
        peak_index_store = np.vstack([peak_index_store, peak_indexes])
        #peak_wavelengths_store = np.vstack([peak_wavelengths_store, peak_wavelengths])
    
        #plt.figure(figsize=(10,6))
        #plt.title(data_set.name[25:])
        #plt.ylabel("Normalised intensity")
        #plt.xlabel("Wavelength(nm)")
        
        #k=0
        #while k < max_peaks:
        #    plt.text(x_red[peak_indexes[k]], y_norm[peak_indexes[k]] + 0.005, 
        #             k)
        #    plt.text(max(x_red)+50, max(y_norm)-((k+1)*0.02), 
        #             "Peak " + str(k) + ": " + str(x_red[peak_indexes[k]]))          
        #    k=k+1
            
       #pplot(x_red, y_norm, peak_indexes_orig)
        
    else:
        pass

#all_info_wl = [wavelength_store, intensity_store, peak_wavelengths_store]
all_info = [wavelength_store, intensity_store, peak_index_store]           
sorting_array = np.concatenate((intensity_store, peak_index_store), axis = 1)
sorted_array = sorted(sorting_array, key=lambda x: x[data_range])

while len(sorted_array[1]) > data_range:
    sorted_array = np.delete(sorted_array, [data_end-data_start], 1)
    
sorted_array = np.delete(sorted_array, 0, 0)

all_sorted_array = np.array(sorted_array, dtype = float)

mapfig = ax2.imshow(all_sorted_array, 
                    extent=[wavelength_start, wavelength_end,
                    0, len(all_sorted_array)], aspect = 'auto', 
                            cmap = parula_map)

del_no = -5
del_front = 5
cut_sorted_array = all_sorted_array

while del_no < del_front:
    cut_sorted_array = np.delete(cut_sorted_array, 0, 0)
    del_no = del_no + 1

mapfig2 = ax3.imshow(cut_sorted_array, 
                    extent=[wavelength_start, wavelength_end,
                    0, len(cut_sorted_array)], aspect = 'auto')
                        
fig2.colorbar(mapfig)
ax2.tick_params(axis='x', labelsize=10)
ax2.tick_params(axis='y', labelsize=10)
ax2.set_xlabel("Wavelength(nm)", fontsize = 10)
ax2.set_ylabel("Spectrum number", fontsize = 10)
ax2.set_title("Normalised DF sorted by peak position", fontsize = 18)

fig3.colorbar(mapfig2)
ax3.tick_params(axis='x', labelsize=10)
ax3.tick_params(axis='y', labelsize=10)
ax3.set_xlabel("Wavelength(nm)", fontsize = 10)
ax3.set_ylabel("Spectrum number", fontsize = 10)
ax3.set_title("Reduced normalised DF sorted by peak position", fontsize = 18)

fig1.savefig(fname[:-3] + "_normalised data.png")
fig2.savefig(fname[:-3] + "_normalised map.png")
fig3.savefig(fname[:-3] + "_reduced map.png")

data_f.close()
