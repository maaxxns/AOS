%How to open a file in Matlab/GNU Octav
%Jean-Michel Ménard 
%24th September 2019

clear all
format long 

readpath=''; %Include the path of the file between the '', for example: readpath='C:\D\PHY7233\', or leave it as '' if the file is in the same folder as this Matlab file
filename1  = 'Interferogram1.txt'; 
    [pos1, I1] = textread([readpath, filename1], '%f %f', 'headerlines',1);


%pos = pos1

%The full stage displacement is 9.166 mm /2, therefore the total
%displacement between the two mirrors of the interferometer is 9.166 mm
pos = pos1./max(pos1)*9.166E-3;
t= pos/3E8;

figure(1)
plot(t,I1,'-b')
title('Interferogram')

%How to do a fourier transform in matlab
%Assuming that the x component has units of time (s):
timestep = mean(diff(t));           % Timestep (s)
Fs = 1./timestep;                   % Sampling frequency in [Hz]    
n_FFT = length(t);
Ef_1 = fft(I1,n_FFT); 

%Keeping only the positive frequencies
NumUniquePts = ceil((n_FFT+1)/2); 
Ef1 = Ef_1(1:NumUniquePts);


%Frequences are defined as:
f = zeros(NumUniquePts,1); %Empty vector
f = (0:NumUniquePts-1)*Fs/n_FFT;    %Defining the frequency vector

figure(2)
plot(f/3E8/100,abs(Ef1),'-b') %The x axis is defined in cm^-1
title('FT(Interferogram)')
axis([0 2*5800 0 1])
xlabel('Wavevector cm^{-1}')
ylabel('Spectral intensity')
    

    
    
    