#!/usr/bin/python
# -*- coding: utf-8 -*-

# © Frank Löffler <knarf@cct.lsu.edu>
# License: GPLv2+

# Compute power spectral density of given data
from plot_defaults import *
from matplotlib.mlab import detrend_linear

xlim = (0,10)
ylim = (-140,-100)

# load data
Fx, Fy = np.loadtxt("hydrobase::rho.maximum.asc.bz2", comments="#", usecols=(1,2), unpack=True)
dt = M_to_ms/(Fx[1]-Fx[0])

# slice data if wanted
start = 0.0
end   = 1.
Fx = Fx[int(start*len(Fx)): int(end*len(Fx))]
Fy = Fy[int(start*len(Fy)): int(end*len(Fy))]

# Plot basics
fig = plt.figure(figsize=(6, 3))
fig.subplots_adjust(top=0.92,bottom=0.16, left=0.11,right=0.98)
ax = fig.add_subplot(1,1,1)

# mode names and frequencies for non-rotating Gamma=2 K=100 TOV stars
modes = {
"F" : 1.44252691528028,
"H1": 3.95408396916149,
"H2": 5.91494894170517,
"H3": 7.77382493405710,
"H4": 9.58768293806276,
"H5": 11.3772129983097,
"H6": 13.1520241905666,
"H7": 14.9172321735655,
}

# plot modes as vertical lines and label them
for mode, freq in modes.items():
  ax.plot((freq,freq), ylim, 'b--')
  ax.text(freq-0.1, ylim[1]+1, mode, fontsize=fontsize)
#ax.plot((2*modes["F"], 2*modes["F"]), (-200,0), 'g--')
#ax.plot((3*modes["F"], 3*modes["F"]), (-200,0), 'g--')
#ax.plot((4*modes["F"], 4*modes["F"]), (-200,0), 'g--')

# plot PSD
frac = 0.5 * len(Fy)
ax.psd(Fy, NFFT=int(frac), pad_to=5*int(frac), noverlap=int(frac*0.5), detrend=detrend_linear, Fs=dt, linestyle='-', color='black', scale_by_freq=False)
#ax.psd(Fy, NFFT=int(len(Fy)*frac), pad_to=10*int(len(Fy)*frac), noverlap=int(len(Fy)*frac*0.50), detrend=detrend_linear, Fs=M_to_ms/(Fx[1]-Fx[0]), linestyle='-', marker='x', color='red')

# plot properties
ax.set_xlim(xlim)
ax.set_ylim(ylim) 

#plt.title('PSD', fontsize=fontsize)
ax.set_xlabel(r'$f$ [kHz]')
ax.set_ylabel(r'PSD [dB/Hz]')

ax.xaxis.set_major_locator(mticker.MaxNLocator(10))
ax.xaxis.set_minor_locator(mticker.MaxNLocator(20))
ax.xaxis.grid(False)
ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
ax.yaxis.set_minor_locator(mticker.MaxNLocator(9))
ax.yaxis.grid(False)
set_tick_sizes(ax, 8, 4)


#plt.show()
plt.savefig('mode_spectrum.pdf')

