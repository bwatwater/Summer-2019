
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import scipy.fftpack
import pandas as pd

'''Code by Bennett Atwater- Summer 2019
Takes csv data and generates a fourier transform. Best for diagnosing noise in long-term signals. '''
#path= file path for data file
path=r'C:\Users\Bennett\Desktop\sensitivity1alt.txt'
data=pd.read_csv(path, sep=" ", header=None, names=['x','y'])
x=data.x
y=data.y
Y = np.fft.fft(y)
freq = np.fft.fftfreq(len(y),.0005)
fig, ax = plt.subplots()
ax.semilogy(freq, np.abs(Y))
pylab.xlim([0,1000])
plt.xticks(np.arange(0, 1000, 20))
plt.xlabel('Frequency (Hz)')
plt.ylabel('Signal Intensity')
plt.show()
