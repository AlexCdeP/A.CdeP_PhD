# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 16:23:38 2018

@author: alexc

"""

from lmfit.models import GaussianModel
import numpy as np
import scipy.constants as const
import peakutils
import matplotlib.pyplot as plt
import scipy.signal as sig
from pandas import DataFrame
from sklearn.cluster import KMeans
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
from matplotlib.collections import LineCollection


class data_manipulators(object):

    def energy_ticks(self, wavelengths):
        E = (((const.h) * 3 * (10 ** (8))) / (wavelengths * (10 ** (-9)))) / (1.6 * 10 ** (-19))
        return ["%.2f" % z for z in E]

    def e_wl_convert(self, energies):
        wls = (((const.h) * 3 * (10 ** (8))) / (energies * (10 ** (-9)))) / (1.6 * 10 ** (-19))
        return wls

    def offset_spectra(self, spectra, s_fact):

        spectra_w_offsets = np.zeros(np.shape(spectra))

        for i in range(len(spectra)):
            spectra_w_offsets[i] = spectra[i] + i * s_fact

        return spectra_w_offsets

    @staticmethod
    def slashfinder(filelocation):
        char_i = 0
        for char in filelocation:
            char_i = char_i + 1
            if char == "/":
                slash_num = char_i
            else:
                pass
        return slash_num

    @staticmethod
    def gaussian_fit(y_data, peak_indexes=None, peak_intensities=None):
        # function assumes that you have numpy imported
        if peak_indexes is None:
            peak_indexes = np.array(np.where(y_data == np.max(y_data)))
        else:
            pass
        if peak_intensities is None:
            peak_intensities = np.array(np.max(y_data))
        else:
            pass

        all_indeces = np.arange(y_data.size)
        fits = []
        for i in range(0, len(peak_indexes)):
            if peak_indexes[i] is not None:
                gauss_centre = peak_indexes[i]
                gauss_max = peak_intensities[i]
                width = np.sqrt(np.abs(np.sum((y_data - gauss_centre) ** 2 * y_data) / np.sum(y_data)))
                fit = lambda t: gauss_max * np.exp(-(t - gauss_centre) ** 2 / (2 * width ** 2))
                fits.append(fit(all_indeces))
            else:
                pass

        return fits

    def gaussian_fit_2(self, y_data, x_data=None, x_peak_indexes=None, y_peak_values=None):

        if x_data is None:
            x_data = np.arange(np.max(y_data))
        else:
            pass
        if x_peak_indexes is None:
            x_peak_indexes = np.array(np.where(y_data == np.max(y_data)))
        else:
            pass
        if y_peak_values is None:
            y_peak_values = np.array([np.max(y_data)])
        else:
            pass

        global multi_gauss_a
        for peak_num in range(len(x_peak_indexes)):
            if peak_num == 0:
                multi_gauss_a = GaussianModel(prefix='g' + str(peak_num) + '_')
            else:
                multi_gauss_a += GaussianModel(prefix='g' + str(peak_num) + '_')

            multi_gauss_a.set_param_hint('g' + str(peak_num) + '_center', value=x_peak_indexes[peak_num],
                                         min=x_peak_indexes[peak_num] - 200, max=x_peak_indexes[peak_num] + 200)

            multi_gauss_a.set_param_hint('g' + str(peak_num) + '_sigma', value=50)

            multi_gauss_a.set_param_hint('g' + str(peak_num) + '_amplitude', value=50.0 * np.sqrt(2 * np.pi) *
                                                                                    y_peak_values[peak_num])

            # FIT ACTUALLY OUTPUTS AREA INSTEAD OF AMPLITUDE

        multi_gauss_fit_a = multi_gauss_a.fit(y_data, x=x_data)  # THIS IS THE LINE ACTUALLY DOING THE FITTING

        return multi_gauss_fit_a

    # The sort coupled modes function sorts data according to peak position. 
    # The peak positions 
    # must be a list of all the peak wavelength values found using a peak finder
    # such as peakutils. The function takes initially deletes first values
    # below a wavelength of 550nm, therefore focussing on peaks at wavelengths 
    # longer than the tranverse mode. It then takes the corresponding intensity
    # values and concantenates them with the peak values, allowing the data
    # to be sorted using the inbuilt sorted function. This sorting also 
    # requires a data range to be specified.

    def sort_coupled_modes(self, spectra, all_ints, all_pk_wls):

        spectra = np.array(spectra)
        l1pks_only = np.zeros(len(spectra))
        for i in range(0, len(all_pk_wls)):
            pkwls_i = all_pk_wls[i]
            pkints_i = all_ints[i]
            l1_arg = np.argwhere(pkints_i == np.max(pkints_i))
            # l1_pkint_i = np.max(pkints_i)
            l1_pkwl_i = pkwls_i[l1_arg]
            l1pks_only[i] = l1_pkwl_i

        sorting_array = np.c_[spectra, l1pks_only]
        sorted_array = np.array(sorted(sorting_array, key=lambda x: x[-1]))
        sorted_array = sorted_array[:, :-1]

        return sorted_array

    def peakfind_DF(self, y_data, x_data, threshold=0.1, peak_width_guess=50,
                    exp_no_pks=3, max_peaks=4, cutoff=0.0, iter=False):

        # Find peak positions using peakutils
        # y_data is the data from the QePro
        # x_data is the wavelengths
        # peak_data_start/end are the range of index values that the search for peaks works on
        # threshold is the fraction of intensity values within y-range which will count as peaks

        # Set peakfinder to only search above set wavelength value with
        # peak_data_start

        # Automatically search for the correct threshold leading the code to
        # find the defined number of peaks

        if iter==False:
            peak_indexes = peakutils.indexes(y_data,
                                             thres=threshold,
                                             min_dist=peak_width_guess * 2)
        else:
            find_peaks = False
            while find_peaks is False and exp_no_pks != 0:
                peak_indexes = peakutils.indexes(y_data,
                                                 thres=threshold,
                                                 min_dist=peak_width_guess * 2)
                if len(peak_indexes) >= exp_no_pks:
                    find_peaks = True
                elif threshold <= 0.01:
                    while not find_peaks:
                        peak_indexes = peakutils.indexes(y_data,
                                                         thres=threshold,
                                                         min_dist=peak_width_guess * 2)
                        if len(peak_indexes) >= exp_no_pks:
                            find_peaks = True
                        elif threshold > 1:
                            print "Didn't find threshold for spectrum"
                            print "Need to change expected number or peaks"
                            break
                        else:
                            threshold += 0.01
                    break
                else:
                    threshold -= 0.01

        peak_indexes = peak_indexes

        if len(peak_indexes) > max_peaks:
            peak_indexes = peak_indexes[:max_peaks]
        else:
            pass

        x_data_low = np.min(x_data) * (1 + cutoff)
        x_data_high = np.max(x_data) * (1 - cutoff)

        peak_wavelengths = []
        peak_intensities = []
        del_arr = []
        for i in range(0, len(peak_indexes)):
            if peak_indexes[i] != 0 and x_data_low < x_data[peak_indexes[i]] < x_data_high:
                peak_wavelengths = np.append(peak_wavelengths,
                                             [x_data[peak_indexes[i]]])
                peak_intensities = np.append(peak_intensities,
                                             [y_data[peak_indexes[i]]])
            else:
                pass

            if peak_wavelengths[0] < 600:
                if peak_intensities[i] < peak_intensities[0]:
                    del_arr.append(i)

        peak_wavelengths = np.delete(peak_wavelengths, del_arr)
        peak_intensities = np.delete(peak_intensities, del_arr)
        peak_indexes = np.delete(peak_indexes, del_arr)

        return peak_indexes, peak_wavelengths, peak_intensities

    def colarr(self, y_data, colours):
        # Check if the y-data array is multi-dimensional
        if np.size(np.shape(y_data)) >= 2:
            normarr = []
            # Make array of maximum normalised values from dataset - assumes peak values within reasonable range 100-900 pnts

            for j in range(0, len(y_data)):
                normj = np.max(y_data[j][100:900]) / np.max(y_data[:, 100:900])
                normarr.append(normj)

            colarr = ((normarr - np.min(normarr)) / np.ptp(normarr)) * (len(colours) - 1)
            colarr = colarr.astype(int)
        # Check if y-data is a vector
        elif np.size(np.shape(y_data)) < 2:
            colarr = y_data / np.max(y_data) * (len(colours) - 1)
            colarr = colarr.astype(int)

        else:
            print 'No array to convert to colour array found'
            colarr = np.linspace(0, 1, len(colours))

        return colarr

    def q_factor_calc(self, x_data, y_data):
        # Find the full width half max for each best_fit
        fwhm_wl_mat = []
        int_ov_fwhm = []
        q_factors = []
        max_intensities = []
        for i in range(len(y_data)):
            # Find the indeces of the two half maxima
            max_intensity = np.max(y_data[i])
            ind_hm_1 = np.min(np.where(y_data[i] > max_intensity / 2))
            ind_hm_2 = np.max(np.where(y_data[i] > max_intensity / 2))
            fwhm_wl = x_data[ind_hm_2] - x_data[ind_hm_1]
            max_freq = const.c / (x_data[np.where(y_data[i] == np.max(y_data[i]))][0] * 10e-9)
            fwhm_freq = const.c / (x_data[ind_hm_1] * 10e-9) - const.c / (x_data[ind_hm_2] * 10e-9)
            fwhm_wl_mat.append(fwhm_wl)
            max_intensities.append(max_intensity)
            int_ov_fwhm.append((np.max(y_data[i]) * 100) / fwhm_wl)
            q_factors.append(max_freq / fwhm_freq)

        return q_factors, max_intensities, fwhm_wl_mat


class plotting_funcs(object):

    # Spectra must be referenced and background first
    def plot_all_df(self, spectra, wls, colours):
        # Will need to use methods from other class in this file
        datman = data_manipulators()
        # set up figure
        fig = plt.figure(1, figsize=(8, 6))
        ax = fig.add_subplot(111)
        # Looks for max within reasonable range i.e. where there are no large numbers due to noise
        maxdf = np.max(spectra[:, 100:900])
        colarr = datman.colarr(y_data=spectra, colours=colours)
        j = 0
        for j in range(0, len(spectra)):
            normj = np.max(spectra[j][100:900]) / np.max(spectra[:, 100:900])
            ax.plot(wls, spectra[j] * 100, alpha=0.15, color=colours[colarr[j]])

        ax.set_xlim(410, 900)
        avdf = np.average(spectra, axis=0)
        ax.set_ylim([0, max(avdf[100:900]) * 100 * 1.5])
        ax.plot(wls, avdf * 100, linewidth=2, color='black')
        ax.set_xlabel("Wavelength (nm)", fontsize=30)
        ax.set_ylabel("Scattering (%)", fontsize=30)
        ax.tick_params(labelsize=25, direction='in', length=10, width=2)

        return fig, ax, avdf

    def plot_all_smooth_df(self, spectra, wls, colours, baseline_subtract=0, start_wl=410, end_wl=900, order=1):
        datman = data_manipulators()

        start_ind = np.where(wls.astype(int) == start_wl)
        start_ind = start_ind[0][0]
        end_ind = np.where(wls.astype(int) == end_wl)
        end_ind = end_ind[0][0]
        wls_filt = wls[start_ind:end_ind]

        fig = plt.figure(1, figsize=(8, 6))
        ax = fig.add_subplot(111)
        filtered_spectra = []

        colarr = datman.colarr(y_data=spectra, colours=colours)
        for j in range(0, len(spectra)):
            filt_data = sig.savgol_filter(spectra[j][start_ind:end_ind], 23, order)
            if baseline_subtract != 0:
                base = peakutils.baseline(filt_data, baseline_subtract)
                ax.plot(wls_filt, (filt_data - base) * 100, alpha=0.15, color=colours[colarr[j]])
                filt_data = filt_data - base
                if np.min(filt_data) < 0:
                    filt_data = filt_data - np.min(filt_data)
                else:
                    pass

                filtered_spectra.append(filt_data)
            else:
                ax.plot(wls_filt, filt_data * 100, alpha=0.15, color=colours[colarr[j]])
                filtered_spectra.append(filt_data)

        avdf_filt = np.average(filtered_spectra, axis=0)
        std_filt = np.std(filtered_spectra, axis=0)
        ax.set_xlim(start_wl, end_wl)
        ax.set_ylim(0, np.max(filtered_spectra[:]) * 100 * 1.2)
        avplusdev = avdf_filt + std_filt
        avmindev = avdf_filt - std_filt
        ax.plot(wls_filt, avdf_filt * 100, linewidth=2, color='black')
        ax.plot(wls_filt, avmindev * 100*0.9, linewidth=2, alpha=0.3, linestyle='dashed', color='red')
        ax.plot(wls_filt, avplusdev * 100*0.9, linewidth=2, alpha=0.3, linestyle='dashed', color='red')
        ax.fill_between(wls_filt, avplusdev * 100, avmindev * 100, color='red', alpha=0.1)
        ax.set_xlabel("Wavelength (nm)", fontsize=30)
        ax.set_ylabel("Scattering (%)", fontsize=30)
        ax.tick_params(labelsize=25, direction='in', length=10, width=2)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.tight_layout(pad=1.6)

        return fig, ax, filtered_spectra, avdf_filt, wls_filt, std_filt

    def plot_all_peaks(self, spectra, wls, colours, threshold=0.1, start_wl=450, end_wl=899, cutoff=0):

        datman = data_manipulators()

        if len(wls) == len(spectra[0]):
            print "spectra same length"
        else:
            print "Oh dear..."

        start_ind = np.where(wls.astype(int) == start_wl)
        start_ind = start_ind[0][0]
        end_ind = np.where(wls.astype(int) == end_wl)
        end_ind = end_ind[0][0]

        # Find peak positions using Alexsdatafuncs peakfinder

        pk_wls_mat = []
        pk_ints_mat = []
        pk_indexes_mat = []
        allpkwls = []
        allpkints = []

        fig = plt.figure(1, figsize=(8, 6))
        ax = fig.add_subplot(111)

        for i in range(len(spectra)):
            pk_inds, pk_wls, pk_ints = datman.peakfind_DF(x_data=wls[start_ind:end_ind],
                                                                                  y_data=spectra[i]
                                                                                  [start_ind:end_ind],
                                                                                  threshold=threshold, cutoff=cutoff,
                                                                                  iter=False)
            pk_wls_mat.append(pk_wls)
            pk_ints_mat.append(pk_ints)
            pk_indexes_mat.append(pk_inds)
            allpkwls = np.append(allpkwls, pk_wls)
            allpkints = np.append(allpkints, pk_ints)

        #Finds the peaks which have maximum intensity for each spectrum and assumes this is the dipole antenna mode

        maxpk_wls = np.zeros(len(pk_wls_mat))
        maxpk_ints = np.zeros(len(pk_wls_mat))
        for i in range(len(pk_wls_mat)):
            arg_max_ind = np.argmax(pk_ints_mat[i])
            if pk_wls_mat[i][arg_max_ind] > 600:
                maxpk_wls[i] = pk_wls_mat[i][arg_max_ind]
                maxpk_ints[i] = pk_ints_mat[i][arg_max_ind]
            else:
                pass

        #Try a different classification method to identify dipole antenna modes which are split or not split

        # for i in range(len(pk_wls_mat)):
        #     if len(pk_wls_mat[i]) <= 2 and pk_wls_mat[i][0] < 600:
        #         l1_only_wls = np.append(l1_only_wls, pk_wls_mat[i], axis=0)
        #         l1_only_ints = np.append(l1_only_ints, pk_ints_mat[i], axis=0)
        #     elif len(pk_wls_mat[i]) > 2 and pk_wls_mat[i][1] < 600:
        #         l1_only_wls = np.append(l1_only_wls, pk_wls_mat[i], axis=0)
        #         l1_only_ints = np.append(l1_only_ints, pk_ints_mat[i], axis=0)
        #     elif len(pk_wls_mat[i]) > 2:
        #         multi_pk_wls = np.append(multi_pk_wls, pk_wls_mat[i], axis=0)
        #         multi_pk_ints = np.append(multi_pk_ints, pk_ints_mat[i], axis=0)
        #     elif len(pk_wls_mat[i]) == 2 and pk_wls_mat[i][0] >= 600:
        #         multi_pk_wls = np.append(multi_pk_wls, pk_wls_mat[i], axis=0)
        #         multi_pk_ints = np.append(multi_pk_ints, pk_ints_mat[i], axis=0)
        #     else:
        #         pass

        ax.scatter(allpkwls, allpkints*100, marker='x', color=colours[80], label='Other pks', alpha=0.7)
        ax.scatter(maxpk_wls, maxpk_ints * 100, marker='o', color=colours[0], label='l=1', alpha=0.6)
        ax.plot(wls, np.average(spectra, axis=0)*100, linestyle="--", color=colours[40])
        ax.set_xlabel("Wavelength (nm)", fontsize=30)
        ax.set_ylabel("Scattering (%)", fontsize=30)
        ax.tick_params(labelsize=25, direction='in', length=10, width=2)
        ax.set_ylim(0, np.max(allpkints)*1.05*100)
        ax.set_xlim(410, 900)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        plt.tight_layout(pad=1.6)

        return fig, ax, pk_wls_mat, pk_ints_mat, allpkwls, allpkints, maxpk_ints, maxpk_wls

    def fit_test(self, spectrum_filt, wls_filt, threshold=0.1, start_wl=450, end_wl=899):
        datman = data_manipulators()

        start_ind = np.where(wls_filt.astype(int) == start_wl)[0][0]
        end_ind = np.where(wls_filt.astype(int) == end_wl)[0][0]

        # Check fit using average filtered data
        peak_indexes, peak_wavelengths, peak_intensities = datman.peakfind_DF(x_data=wls_filt[start_ind:end_ind],
                                                                              y_data=spectrum_filt[start_ind:end_ind],
                                                                              threshold=threshold)

        fig = plt.figure(1, figsize=(7, 7))
        ax = fig.add_subplot(111)
        ax.scatter(peak_wavelengths, peak_intensities * 100)
        ax.plot(wls_filt, spectrum_filt * 100)
        ax.set_ylim(0, max(spectrum_filt) * 100 * 1.05)
        ax.set_xlim(410, 900)
        ax.tick_params(labelsize=22, direction='in', length=10, width=2)
        ax.set_ylabel("Scattering (%)", fontsize=25)
        ax.set_xlabel("Wavelength (nm)", fontsize=25)

        return fig, ax

    def peak_histogram(self, allpkwls, wls, colours, av_spectrum=None, ra=[420, 900], bins=50):
        datman = data_manipulators()
        if av_spectrum is None:
            av_spectrum = np.zeros(len(wls))
        else:
            pass

        # make histogram of peak positions and intensities
        fig = plt.figure(1, figsize=(8, 6))
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()
        ax2.plot(wls, av_spectrum * 100, color='black', linestyle='--')
        n, bins, patches = ax.hist(allpkwls, range=ra, bins=bins, alpha=0.5)
        colarr = datman.colarr(n, colours)

        patchnum = 0
        for patch in patches:
            patch.set_facecolor(color=colours[colarr[patchnum]])
            patchnum = patchnum + 1

        ax.set_xlim(400, 900)
        ax.set_ylim(0, 70)
        ax.set_ylabel('Peak count', fontsize=30)
        ax.tick_params(labelsize=25)
        ax.set_xlabel('Wavelength (nm)', fontsize=30)
        ax2.set_ylabel('Scattering (%)', fontsize=30, rotation=-90, labelpad=28)
        ax.tick_params(labelsize=22, direction='in', length=10, width=2)
        ax2.tick_params(labelsize=22, direction='in', length=10, width=2)
        ax2.set_ylim(0.001, max(av_spectrum) * 100 * 1.05)
        plt.tight_layout(pad=1.6)

        return fig, ax, ax2, n, bins

    def kmeans_pkplot(self, pkwls_vec, pkints_vec, n_clusters=3):
        # k means clustering

        data = {'x': pkwls_vec, 'y': pkints_vec}
        dframe = DataFrame(data, columns=['x', 'y'])

        kmeans = KMeans(n_clusters=n_clusters).fit(dframe)
        centroids = kmeans.cluster_centers_
        labels = kmeans.labels_

        fig = plt.figure(1, figsize=(8, 6))
        ax = fig.add_subplot(111)
        ax.scatter(dframe['x'], dframe['y'] * 100, c=kmeans.labels_.astype(float), s=50, alpha=0.5)
        ax.scatter(centroids[:, 0], centroids[:, 1] * 100, c='red', s=50)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


        return fig, ax, dframe, kmeans, centroids, labels

    def default_hist(self, data_vector, colours='black', range=None, bins=None):
        datman = data_manipulators()
        if range is None:
            range = [np.min(data_vector), np.max(data_vector)]
        else:
            pass
        if bins is None:
            bins = (len(data_vector) / max(data_vector)).astype(int) * 2
        else:
            pass

        fig = plt.figure(1, figsize=(7, 7))
        ax = fig.add_subplot(111)
        n, bins, patches = ax.hist(data_vector, range=range, bins=bins)
        colarr = datman.colarr(n, colours)

        patchnum = 0
        for patch in patches:
            patch.set_facecolor(color=colours[colarr[patchnum]])
            patchnum = patchnum + 1

        ax.set_ylabel("Counts", fontsize=25)
        ax.set_xlabel("vector_name", fontsize=25)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_xlim(np.min(data_vector), np.max(data_vector))

        return fig, ax, n, bins

    def offset_map(self, spectra, wls, lthickness=40, offset_fact=0.025, norm=True):

        datman = data_manipulators()


        offset = np.max(spectra) * offset_fact
        offsetted_spectra = datman.offset_spectra(spectra, s_fact=offset)

        fig1 = plt.figure(1, figsize=(6, 11))
        ax1 = fig1.add_subplot(111)

        colours = cm.viridis(np.linspace(0, 1, len(offsetted_spectra)))
        for i in range(0, len(offsetted_spectra) - 1):
            lwidths = (offsetted_spectra[i] - i * offset) * lthickness
            points = np.array([wls, offsetted_spectra[i]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc_2 = LineCollection(segments, linewidths=lwidths, color='black')
            ax1.add_collection(lc_2)
            spec1 = np.array(offsetted_spectra[i], dtype=float)
            spec2 = np.ones(len(wls)) * i * offset
            ax1.fill_between(wls, spec1, spec2, color=colours[i],
                             alpha=0.1, zorder=3)

        ax1.set_xlim(450, np.max(wls))
        ax1.set_xlabel("Wavelength (nm)", fontsize=30)
        ax1.tick_params(axis='x',labelsize=25, direction='in', length=10, width=2)

        ax1.set_ylim(0, )
        ax1.set_yticklabels([])

        return fig1, ax1, offsetted_spectra

# class test():

#    def __init__(self, model, colour):
#        self.model = model
#        self.colour = colour

#    def default(self, spectra):
#        datman = data_manipulators()

#        dataman.colarr(spectra) = self.colour
#        self.ylim =
