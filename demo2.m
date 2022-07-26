% ==================================================
% Demo Script (2)
%
% Objectives:
% --------------------------------------------------
% Use the Signal MultiResolution Analyzer App
% and the Wavelet Analyzer to decompose and then
% reconstruct an EEG recording.

clear all; clc;

% ==================================================
% Read EEG recordings from the EDF file
% ==================================================

X = edfread("SN001.edf");      % Read the first EDF file

t0 = duration("00:00:00");     
dt = duration("00:00:30");
range = timerange(t0,t0+dt);

sig = X{range,"EEGF4_M1"};     % extract the first 30 seconds
                               % from the F4-M1 EEG channel

sig = cell2mat(sig);           % cell-array => 1D array

% ==================================================
% 2) Perform MRA Decomposition and Reconstruction
% ==================================================

% Logical array for selecting reconstruction elements
levelForReconstruction = [
    false, ...     % discard details at level 1    (scale: 64Hz-128Hz)
    false, ...     % discard details at level 2    (scale: 32Hz-64Hz)
    false, ...     % discard details at level 3    (scale: 16Hz-32Hz)
    false, ...     % discard details at level 4    (scale: 8Hz-16Hz)
    true,  ...     % keep details at level 5       (scale: 4Hz-8Hz)
    true,  ...     % keep details at level 6       (scale: 2Hz-4Hz)
    true,  ...     % keep details at level 7       (scale: 1Hz-2Hz)
    true];         % keep approximation at level 7 (scale: 0Hz-1Hz)

% Perform the decomposition using modwt
wt = modwt(sig,'db5',7);

% Construct MRA matrix using modwtmra
mra = modwtmra(wt,'db5');

% Sum down the rows of the selected multiresolution signals
sig1 = sum(mra(levelForReconstruction,:),1);

% ==========================================================
% 3) Plot the results
% ==========================================================

% Add a time axis
t0 = 0; dt = 30; fs = 256;
time = linspace(t0,t0+dt,dt*fs);

% Original (sig) vs Reconstructed (sig1)
figure(1); 
plt = plot(time,sig, time,sig1);
plt(1).LineWidth = 0.5;
plt(2).LineWidth = 2;
xlabel('time in seconds');
ylabel('Amplitude in microvolts');
title('Original vs Reconstructed Signal');
legend("original", "reconstructed");

% =========================================================
% 4) 2nd, 3rd and 4th order statistics on all frequency
% scales
% =========================================================

var = std(mra,0,2);       % second order statistics on all scales
skw = skewness(mra,0,2);  % third order statistics on all scales
krt = kurtosis(mra,0,2);  % fourth order statistics on all scales

scales = [128, 64, 32, 16, 8, 4, 2, 1, 0];
num_of_scales = 8;

for i = 1:1:num_of_scales
    f1 = scales(i+1); f2 = scales(i);
	
    fprintf("Frequency scale: %dHz-%dHz\n",f1,f2);
    fprintf("variance = %f\n", var(i));
    fprintf("skewness = %f\n", skw(i));
    fprintf("kurtosis = %f\n", krt(i));
    fprintf("==========================\n");
    fprintf("\n");
end
