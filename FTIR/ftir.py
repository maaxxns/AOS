import numpy as np 
import matplotlib.pyplot as plt 
from scipy.fft import fft, fftfreq
from scipy.constants import c
#functions
def time_step(pos):
    pos = pos/np.max(pos) * 9.166*10**(-3)
    return pos/(3*10**8)

def FFT_func(I, t):
    N = len(t) #length of t1
    timestep = t[2]-t[3]
    FX = fft(I) #the fourier transform of the intensity. 
    FDelay = fftfreq(N, d=timestep) #FFT of the time to frequencies. 
    return [FDelay, FX]


#First get the data from txt
data1  = np.genfromtxt('Interferogram1.txt',skip_header=1, delimiter='\t')
data2 = np.genfromtxt('Interferogram2.txt',skip_header=1, delimiter='\t')
data_ref =  np.genfromtxt('Interferogram_Ref.txt',skip_header=1, delimiter='\t')

#Lets put the data into arrays for position and intensity
pos1 = data1[:,0]
I1 = data1[:,1]
pos2 = data2[:,0]
I2 = data2[:,1]
pos_ref = data_ref[:,0]
I_ref = data_ref[:,1]

#Now we calculate the time difference per step
t1 = time_step(pos1)
t2 = time_step(pos2)
t_ref = time_step(pos_ref)

#plot the intensity against time delay
plt.figure()
plt.plot(t1*10**(12), I1, label='Intensity')
plt.xlabel(r'$\delta t / ps$')
plt.ylabel('Intensity / W')
plt.grid()
plt.savefig('interferogram1.pdf')
plt.close()
########################################################################################################################

#now fourier transform with the function that is defined at the top

freq1, amp1 = FFT_func(I1, t1)
freq_, amp_ref = FFT_func(I_ref,t_ref)

plt.figure()
plt.plot(np.abs(freq1*10**(-12)), np.abs(amp1),linewidth=0.1, label='FFT') #now we plot. because the FX is imaginary we take its absolute value.
plt.xlabel(r'$f / THz$')
plt.ylabel('Spectral Intensity')
plt.legend()
plt.savefig('F[interferogram1].pdf')
plt.close()
#########################################################################################################################

I_trans = (I1)

#t1 and t_ref are the same
#wavevector k = w/c nu_bar = nu/c

freq_trans, amp_trans = FFT_func(I_trans, t_ref) #calculate FFT for the transmission

amp_trans = amp_trans/amp_ref

nu_bar = np.abs(freq_trans/(c))/100 #calculate the wavenumber

index = np.where((nu_bar>=500) & (nu_bar<5800))[0] #search for the wavenumber 500<v<5800 and use those as index

plt.figure()
plt.plot(np.abs(nu_bar[index]), np.abs(amp_trans[index]),linewidth=0.1, label='FFT') #now we plot. because the FX is imaginary we take its absolute value.
plt.xlabel(r'$\bar{\nu} / cm$')
plt.ylabel('Spectral Tranmission Intensity')
plt.legend()
plt.savefig('F[interferogram_trans].pdf')

#########################################################################################################################
#average in time domain


I_trans_avg = ((I1+I2)/2) #average with both datasets and than calculate transmission

freq_trans_avg, amp_trans = FFT_func(I_trans_avg, t_ref) #calculate FFT for the transmission

amp_trans = amp_trans/amp_ref

nu_bar_avg = np.abs(freq_trans_avg/(c))/100 #calculate the wavenumber

index2 = np.where((nu_bar_avg>=500) & (nu_bar_avg<5800))[0] #search for the wavenumber 500<v<5800 and use those as index

plt.figure()
plt.plot(np.abs(nu_bar_avg[index2]), np.abs(amp_trans[index2]),linewidth=0.1, label='FFT') #now we plot. because the FX is imaginary we take its absolute value.
plt.xlabel(r'$\bar{\nu} / cm$')
plt.ylabel('Spectral Tranmission Intensity')
plt.legend()
plt.savefig('F[interferogram_trans_avg_in_time].pdf')

##########################################################################################################################
#average in spectral domain

freq_ref, amp_ref = FFT_func(I_ref, t_ref) #calculate FFT for the transmission

f1, amp_1 = FFT_func(I1, t1) #calculate FFT for the transmission

f2, amp_2 = FFT_func(I2, t2) #calculate FFT for the transmission

amp_trans_avg = (amp_1 + amp_2)/2

amp_trans_avg = amp_trans_avg/amp_ref


nu_bar_avg = np.abs(freq_trans_avg/(c))/100 #calculate the wavenumber

index2 = np.where((nu_bar_avg>=500) & (nu_bar_avg<5800))[0] #search for the wavenumber 500<v<5800 and use those as index

plt.figure()
plt.plot(np.abs(nu_bar_avg[index2]), np.abs(amp_trans_avg[index2]),linewidth=0.1, label='FFT') #now we plot. because the FX is imaginary we take its absolute value.
plt.xlabel(r'$\bar{\nu} / cm$')
plt.ylabel('Spectral Tranmission Intensity')
plt.legend()
plt.savefig('F[interferogram_trans_avg_in_freq].pdf')


#############################################################################################################################
#spectral absorption 


amp_abs = 1-np.abs(amp_trans_avg) 
lines = np.array([550, 750, 1630, 1650, 3200, 3650])
plt.figure()
plt.plot(np.abs(nu_bar_avg[index2]), np.abs(amp_abs[index2]),linewidth=0.1, label='FFT') #now we plot. because the FX is imaginary we take its absolute value.
for i in range(len(lines)):
    plt.vlines(lines[i], 0.0, 1, colors='r', linestyle='dashed', linewidth=0.2)
plt.xlabel(r'$\bar{\nu} / cm$')
plt.ylabel('Spectral Absorption Intensity')
plt.legend()
plt.savefig('F[interferogram_abs_avg_in_freq].pdf')
