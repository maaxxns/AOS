%Help file for CS#3
%Jean-Michel Ménard 
%12th January 2023

clear all
format long 

cd('C:\Users\'); %Replace by your folder containing the 2 xlsx files
reference = xlsread('Unknown_Ref');
sample = xlsread('Unknown_Sample');

thickness = 0.9E-3; %[m]
        
%REFERENCE
yr(1,:) = reference(:,2);
xr(1,:) = reference(:,1);


%SAMPLE
ys(1,:) = sample(:,2);
xs(1,:) = sample(:,1);


%% 
% You can use this section to plot the reference and sample time-domain data in the same graph.


%%
% You can use this section to apply the time-domain filter to remove "artifacts" in the time-domain data and plot the resulting traces in the same figure. (make sure that after filtering, you have the same number of data points
% for each set of data, so that you will have the same number of data points after performing the Fourier Transform)


%% 
% You can use this section to perform the Fourier Transform and plot the specta of both the reference and sample measurements in the same figure.

% Example of Fourier Tansfrom applied to reference function
% Define frequency array
t = reference(:,1); % [ps]
ref = reference(:,2);
n = length(t);
dt = mean(diff(t)); %This is just one way of obtaining this value, one could use a different calculation.
dV = 1/dt; 
V = (0:n-1)*dV/n; % [THz] This is the vector for the Frequency

% Perform Fourier Transform 
refFT = fft(ref);
%Keeping only the positive frequencies
positivePts = ceil((n+1)/2); 
ref = ref(1:positivePts);

% Plot the spectrum of the reference measurements
figure(3) 
plot(V,abs(refFT))
xlabel('Frequency [THz]')
ylabel('Electric field [arb]')
xlim([0 2])
% Perfrom Fourier Transform: 

%%
% You can use this section to calculate the refractive index and absorption coefficient of the material from 0.2-0.9 THz and plot them.

% Below is an example of extracting the phase from the FFT function. Two matlab function could be used here: 
% a. angle (theta = angle(z) returns the phase angle in the interval [-pi,pi] for each element of a complex array z.)
% b. unwrap (Whenever the jump between consecutive angles is greater than or equal to pi radians, unwrap recalculate the angles 
%    by adding multiples of ±2pi until the jump is less than pi = 3.1415926...)

H0 = spectrum_sample./spectrum_reference;
% === Extract phase info
phase = unwrap(angle(H0));


