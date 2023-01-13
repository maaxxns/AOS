import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.collections import LineCollection
import matplotlib.colors as clr
from functions import *
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap

phi_min = -50 #Minimum rotation angle of the grating: in deg
phi_max = 50 #Maximum rotation angle of the grating: in deg
phi_acc = 10000 #Accuracy of the vector phi
phi = np.linspace(phi_min,phi_max,phi_acc)*np.pi/180 #Rotation of the graing in [rad] 

#Write your code here

#known variables


d = 3.7*10**(-6) #m
N = 500
sigma = 20*np.pi/180# rad
beta = (30)*np.pi/180-phi #rad
alpha = -phi
m=1
f = 0.3

lambda_min = 750*10**(-9) #Minimum wavelength of the incoming beam (microns)
lambda_max = 1100*10**(-9) #Maximum wavelength of the incoming beam (microns)
lambda_acc = 100 #Accuracy of the vector wavelength
lambda_ = np.linspace(lambda_min,lambda_max,lambda_acc) #

#plotting

fig = plt.figure(figsize=(5,20))

plt.xlabel(r'$\phi \,/\, °$')
plt.ylabel('Diffraction effiency / %')

#ax_2.set_xlabel(r'$\phi \,/\, °$')
#ax_2.set_ylabel('Diffraction effiency / %')

norm = Normalize(vmin=lambda_min, vmax=lambda_max)
cmap = get_cmap('jet')
print('alpha, ', np.arcsin(m*lambda_min/d - np.sin(beta))*180/np.pi)
#print(slitsize(d, 800*10**(-9), alpha, f))
Trans = np.zeros(phi_acc)
for i in range(len(lambda_)):
    Trans =  T_fct(lambda_[i],  N, d, alpha, beta, sigma) + Trans

    plt.plot(phi*180/np.pi, T_fct(lambda_[i],  N, d, alpha, beta, sigma),c = cmap(norm(lambda_[i])))
    #ax_2.plot(phi*180/np.pi, T_fct(lambda_[i],  N, d/2, alpha, beta, sigma),c = cmap(norm(lambda_[i])))

plt.title('d='+ str(d) + r'$\mu m$' + ', N=' + str(N) )
#ax_2.set_title('d=0.75' + r'$\mu m$')
#plt.title('750nm to 1100nm, d=1.5' + r'$\mu m$' + r'$\sigma = 20°, \alpha = 0 °, \beta = 30°$')
plt.tight_layout()
plt.ylim(0,1)
'''plt.legend()
ax_2.legend()'''
plt.savefig('spectrometer.pdf')
plt.show()


#print('max slit size: ', np.array(slitsize(d, 800*10**(-9), alpha, f)).max(), ' min slitsize: ', np.array(slitsize(d, 800*10**(-9), alpha, f)).min())
plt.close()
plt.figure()
plt.xlabel(r'$\alpha \,/\, °$')
plt.ylabel('slitsize / m')
plt.yscale('log')
plt.plot(alpha*180/np.pi, slitsize(d, 800*10**(-9), alpha, f), label='slitsize')
plt.grid()
plt.legend()
plt.savefig('slitsize.pdf')