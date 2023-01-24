import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import uncertainties as un
from matplotlib.cm import get_cmap
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
from scipy.constants import c

#Read the excel file
data_ref = pd.read_excel('Unknown_Ref.xlsx')
data_sam = pd.read_excel('Unknown_Sample.xlsx')
data_ref = np.array(data_ref)
data_sam = np.array(data_sam)
data_ref[:,0] = data_ref[:,0]*10**(-12)
data_sam[:,0] = data_sam[:,0]*10**(-12)

#plot the intensity against time delay
plt.figure()
plt.plot(data_ref[:,0]*10**(12), data_ref[:,1], label='Reference')
plt.plot(data_sam[:,0]*10**(12), data_sam[:,1], label='Sample')
plt.xlabel(r'$ t/ps $')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.title('The reference and sample data set')
plt.savefig('THz1.pdf')
plt.close()

#############################################################

def super_gaussian(x, amplitude=1.0, center=0.0, sigma=1, expon=1):
    """super-Gaussian distribution
    super_gaussian(x, amplitude, center, sigma, expon) =
        (amplitude/(sqrt(2*pi)*sigma)) * exp(-abs(x-center)**expon / (2*sigma**expon))
    """
    sigma = max(1.e-15, sigma)
    return ((amplitude/(np.sqrt(2*np.pi)*sigma))
            * np.exp(-abs(x-center)**expon / 2*sigma**expon))

filter_ref = super_gaussian(x=data_ref[:,0]*10**(12), center=10**(12)*data_ref[find_peaks(data_ref[:,1], height=0.2)[0],0])
filter_sam = super_gaussian(x=data_sam[:,0]*10**(12), center=10**(12)*data_sam[find_peaks(data_sam[:,1], height=0.1)[0],0])

print('determined peak of the ref data: ', data_ref[find_peaks(data_ref[:,1], height=0.2)[0],0])
print('determined peak of the sam data: ', data_sam[find_peaks(data_sam[:,1], height=0.1)[0],0])

plt.figure()
#plt.plot(data_ref[:,0], data_ref[:,1], label='Reference')
plt.plot(data_ref[:,0]*10**(12), data_ref[:,1]*filter_ref, label='Reference Filtered')
plt.plot(data_sam[:,0]*10**(12), data_sam[:,1]*filter_sam, label='Sample Filtered')
plt.plot(data_ref[:,0]*10**(12), filter_ref, label='super Gaussian for Reference')
plt.plot(data_ref[:,0]*10**(12), filter_sam, label='super Gaussian for sample')
plt.xlabel(r'$ t/ps $')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.title('The filtered data sets')
plt.savefig('THz3_gauss.pdf')
plt.close()

#x_0_ref = np.argwhere(data_ref == 0)[0]
#x_0_sam = np.argwhere(data_sam >=  7.2)[0]

plt.figure()
#plt.plot(data_ref[:,0], data_ref[:,1], label='Reference')
plt.plot(data_ref[43:153, 0]*10**(-12), data_ref[43:153,1], label='Reference Filtered')
plt.plot(data_sam[80:190, 0]*10**(-12), data_sam[80:190,1], label='Sample Filtered')
plt.xlabel(r'$ t/ps $')
plt.ylabel('Amplitude')
plt.legend()
plt.grid()
plt.title('The filtered data sets')
plt.savefig('THz3_1.pdf')
plt.close()

data_ref_gaus = data_ref[:,1]*filter_ref
ref_time = data_ref[:, 0]
data_sam_gaus = data_sam[:,1]*filter_sam
sam_time = data_sam[:, 0]

data_ref = data_ref[43:153] #apply filter
data_sam = data_sam[80:190]

#########################################################################################
def FFT_func(I, t): 
    N = len(t) #length of t1
    timestep = np.abs(t[2]-t[3])
    FX = fft(I)[:N//2] #the fourier transform of the intensity. 
    FDelay = fftfreq(N, d=timestep)[:N//2] #FFT of the time to frequencies. 
    return [FDelay, FX]

freq_ref_gaus, amp_ref_gaus = FFT_func(data_ref_gaus, ref_time)
freq_sam_gaus, amp_sam_gaus = FFT_func(data_sam_gaus, sam_time)

freq_ref, amp_ref = FFT_func(data_ref[:,1], data_ref[:,0])  #in Hz
freq_sam, amp_sam = FFT_func(data_sam[:,1], data_sam[:,0])

plt.figure()
plt.plot(freq_ref*10**(-12), np.abs(amp_ref), label='Reference filtered FFT')
plt.plot(freq_sam*10**(-12), np.abs(amp_sam), label='Sample filtered FFT')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel('Spectral Amplitude')

plt.legend()
plt.grid()
plt.title('The FFT of the filtered data sets')
plt.savefig('THz4_1.pdf')
plt.close()

plt.figure()
plt.plot(freq_ref_gaus*10**(-12), np.abs(amp_ref_gaus), label='Reference filtered FFT')
plt.plot(freq_sam_gaus*10**(-12), np.abs(amp_sam_gaus), label='Sample filtered FFT')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel('Spectral Amplitude')
plt.xlim(0, 2)
plt.legend()
plt.grid()
plt.title('The FFT of the filtered data sets')
plt.savefig('THz4_gauss.pdf')
plt.close()

##############################################################################################

freq_ref_gaus = freq_ref_gaus[1:] # we need to have same shapes
H_0 = amp_sam_gaus/amp_ref_gaus
angle = np.angle(H_0[1:]) #angle between complex numbers
phase = (np.unwrap(angle,period=2*np.pi))  #phase 
d = 0.9*10**(-3) #millimeter

n = 1 - c/(freq_ref_gaus* d) *phase
n_im_gauss = c/(freq_ref_gaus *d) *(np.log((4*n)/(n-1)**2)) - np.log(np.abs(phase))

plt.figure()
#plt.plot(freq_ref_gaus*10**(-12),phase, label='phase')
plt.plot(freq_ref_gaus*10**(-12),n_im_gauss, label='imagenary part of refractive index')
plt.plot(freq_ref_gaus*10**(-12), n, label='real part of refractive index')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel('n (arb.)')
plt.xlim(0.18, 1)
plt.legend()
plt.grid()
plt.title('The real part of the refractive index')
plt.savefig('THz5_gaus.pdf')
plt.close()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_ref = freq_ref[1:] # we need to have same shapes
H_0 = amp_sam/amp_ref

angle = np.angle(H_0[1:]) #angle between complex numbers
phase = (np.unwrap(angle,period=2*np.pi))  #phase 
d = 0.9*10**(-3) #millimeter

n = 1 - c/(freq_ref* d) *phase
n_im = c/(freq_ref *d) *(np.log((4*n)/(n-1)**2)) - np.log(np.abs(phase))

plt.figure()
#plt.plot(freq_ref*10**(-12),phase, label='phase')
plt.plot(freq_ref*10**(-12),n_im, label='imagenary part of refractive index')
plt.plot(freq_ref*10**(-12), n, label='real part of refractive index')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel('n (arb.)')
plt.xlim(0.18, 1)
plt.legend()
plt.grid()
plt.title('The real part of the refractive index')
plt.savefig('THz5_1.pdf')
plt.close()


##################################################################################################
alpha_gauss = 2*freq_ref_gaus *n_im_gauss/c 

plt.figure()
plt.plot(freq_ref_gaus*10**(-12), alpha_gauss*10**(-2), label='Absorption coefficient')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel(r'$\alpha$')
plt.xlim(0.18, 1)
plt.legend()
plt.grid()
plt.title('The absorption coeffiecient')
plt.savefig('THz6_gauss.pdf')
plt.close()

alpha = 2*freq_ref *n_im/c 

plt.figure()
plt.plot(freq_ref*10**(-12), alpha*10**(-2), label='Absorption coefficient')
#plt.plot(data_ref[:,0], filter_ref)
plt.xlabel(r'$ \omega/THz $')
plt.ylabel(r'$\alpha$')
plt.xlim(0.18, 1)
plt.legend()
plt.grid()
plt.title('The absorption coeffiecient')
plt.savefig('THz6.pdf')
plt.close()
