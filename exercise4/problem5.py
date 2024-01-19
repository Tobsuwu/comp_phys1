""" This programs calculates the fourier transform from the given signal. """

import numpy as np
import matplotlib.pyplot as plt


# ADD CODE: read in the data here 

data = np.loadtxt('signal_data.txt')
t = data[:,0]
f = data[:,1]

dt = t[1]-t[0]
N = len(t)

# Fourier coefficients from numpy fft normalized by multiplication of dt
F = np.fft.fft(f)*dt
F_filt = np.copy(F)

# frequencies from numpy fftfreq
freq = np.fft.fftfreq(len(F),d=dt)
freq_filt = np.zeros_like(freq)

for i in range(len(freq)):  # Filters the signal so that only frequencies between 40 Hz and 60 Hz remain.
    if 40 < freq[i] < 60:
        freq_filt[i] = freq[i]
    else:
        F_filt[i] = 0

# inverse Fourier with numpy ifft (normalization removed with division by dt).
# We do the same for the filtered signal.

iF = np.fft.ifft(F/dt)
iF_filt = np.fft.ifft(F_filt/dt)

# positive frequencies are given as
# freq[:N//2] from above or freq = np.linspace(0, 1.0/dt/2, N//2)

fig, ax = plt.subplots()
# plot over positive frequencies the Fourier transform
# This figure shows which frequencies (y-axis) have noticeable signal in them. In this case it is frequencies
# 50 Hz, 80 Hz and 120 Hz

ax.plot(freq[:N//2], np.abs(F[:N//2]))
ax.set_xlabel(r'$f$ (Hz)')
ax.set_ylabel(r'$F(\omega/2\pi)$')
 
# plot the "signal" and test the inverse transform
fig, ax = plt.subplots()
ax.plot(t, f,t,iF.real,'r--')
ax.set_xlabel(r'$t$ (s)')
ax.set_ylabel(r'$f(t)$')

# plot the filtered signal with the original signal. The figure seems reasonable. Filtered signal has the same
# shape as original signal but with a cut-off due to filtering.

fig, ax = plt.subplots()
ax.plot(t, f, t, iF_filt.real, 'r--')
ax.set_xlabel(r'$t$ (s)')
ax.set_ylabel(r'$f(t)$')
ax.legend(["Original data", "Filtered data"])

plt.show()
