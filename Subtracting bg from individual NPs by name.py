# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 14:21:34 2018

@author: alexc
"""

import matplotlib.pyplot as plt
 
fig11 = plt.figure(11)
ax11 = fig11. add_subplot(111)
spectra_store = []
for i in range(0,len(all_info_wl[0])-1):
    if all_info_wl[0][i][0][:11] == all_info_wl[0][i+1][0][:11]:
        subtracted_spec = (all_info_wl[2][i] - all_info_wl[2][i+1])
        spectra_store.append(subtracted_spec)
        ax11.plot(x_red, subtracted_spec)
    else:
        pass

fig13 = plt.figure(13)
ax13 = fig13.add_subplot(111)
spectra_store = np.array(spectra_store, dtype = float)
ax13.imshow(spectra_store,
            extent=[wavelength_start, wavelength_end,
            0, len(spectra_store)], aspect = 'auto', 
                    cmap = parula_map,interpolation = 'None')
    
