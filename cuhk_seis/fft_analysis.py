#/usr/bin/env python
# Author: Jinping ZI
# Histroy:
#		2021-03-14 Define plot_fft function
#

import numpy as np
import matplotlib.pyplot as plt

def plot_fft(freq,amp,ylabel=None):
	'''
	Make three plots of amplitude:
		1.Amplitude - Frequency
		2.Amplitude - Log-frequency
		3.Amplitude - Period.
	No zero frequency allowed. Error in the third plot
	'''
	fig1 = plt.figure(1,figsize=(10,4))
	plt.rcParams['xtick.bottom'] = True
	plt.rcParams['xtick.labelbottom'] = True
	plt.rcParams['xtick.top'] = False
	plt.rcParams['xtick.labeltop'] = False

	ax1 = plt.subplot(1,3,1)
	plt.xlabel("Frequency (Hz)")
	plt.ylabel("Amplitude spectrum")
	plt.plot(freq,amp)

	ax2 = plt.subplot(1,3,2)
	plt.xlabel("Frequency (Hz)")
	plt.ylabel("Amplitude spectrum")
	plt.plot(freq,amp)
	plt.xscale("log")

	ax3 = plt.subplot(1,3,3)
	plt.xlabel("Period (s)")
	plt.ylabel("Amplitude spectrum")
	plt.plot(np.reciprocal(freq),amp)

	plt.tight_layout()
	plt.show()

