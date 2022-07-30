% ======================================================
% Demo Script (3)
% ------------------------------------------------------
%
% Author: Christodoulos Michaelides, July 2022
% ------------------------------------------------------
%
% Objectives:
% ------------------------------------------------------
% 1) Use Multi-Resolution-Analysis on the entire 
% recording of an EEG/EOG/ECG channel
%
% 2) Choose the most appropriate frequency scales
% to reconstruct an approximation of the original
% recording
% 
% 3) Estimate 2nd, 3rd and 4th order moments
% for every frequency scale independently. 
% Use a small time window of 1-3 seconds, when 
% estimating variance, skewness and kurtosis in 
% order to obtain localised estimates of the above
% quantities.
% 
% 4) See if there is any correlation between 
% variance/kurtosis/skewness and transitional phenomena
% (K-complex, arousal periods) on the reconstructed EEG
% by visual inspection of the waveforms.
% ======================================================

% reset your workspace and clear terminal
clear all; close all; clc;

% ======================================================
% 1) Parameters of MRA Decomposition and Reconstruction
% (Choose any value that you want to experiment with)
% ======================================================

input_file = "SN001.edf";    % name of input file
channel = 1;                 % index of selected channel
w = duration("00:00:03");    % window length when estimating moments
wavelet = "db5";             % wavelet for MRA 
levels = 7;                  % number of levels on MRA binary tree 
pivot = 5;                   % Use every frequency scale below pivot  
                             % (including the pivot scale) for 
                             % signal reconstruction

% ------------------------------------------------------
% Do not change anything below that point, 
% unless you know exactly what you are doing.
% ------------------------------------------------------

% number of frequency scales
num_of_scales = levels + 1;                

% ======================================================
% 2) Read data from the EDF file
% ======================================================

fprintf("Loading input file ...  ");

% Read the entire EDF file and
% choose the appropriate channel
X = edfread(input_file);
sig = X{:,channel};    

% Read file metadata
info = edfinfo(input_file);

% N:  number of data records per recording
% d:  duration of every data record in seconds
% n:  samples per data record
% fs: sampling frequency in Hertz
N = info.NumDataRecords;
d = seconds(info.DataRecordDuration);
n = info.NumSamples(channel);
fs = n / d;

% construct a time axis for
% the EEG/EOG/ECG channel
time = linspace(0, N*d, N*d*fs);

% cell-array -> 1D array
sig = cell2mat(sig);       

% delete unused variables
clear X;

fprintf("Done\n\n");

% ======================================================
% 3) Perform MRA decomposition and reconstruction
% ======================================================

fprintf("Performing MRA ... ");

% Logical array for selecting reconstruction scales
levelForReconstruction = (1:1:num_of_scales) >= pivot;

% Perform the decomposition using modwt
wt = modwt(sig,wavelet,levels);

% Construct MRA matrix using modwtmra
mra = modwtmra(wt,wavelet);

% Sum down the rows of the selected multiresolution signals
sig1 = sum(mra(levelForReconstruction,:),1);

% Delete unused variables
clear wt;  

fprintf("Done\n\n");

% ======================================================
% 4) Plot of Original vs Reconstructed Waveform
% ======================================================            

% Reconstructed vs Original Signal
figure(1); 

plt = plot(time,sig, time,sig1);
plt(1).LineWidth = 0.5; 
plt(2).LineWidth = 2.0;

xlabel("time in seconds");
ylabel("Amplitude in microVolts");
title("Original vs Reconstruction");

% ======================================================
% Estimate variance/skewness/kurtosis
% in a sliding window.
% ======================================================

% length of sliding window in samples
w = seconds(w);
l = w * fs;

% number of segments created by sliding window
K = N * d * fs / l;

% Plots of selected frequency scales
% and sliding window statistics
for i = pivot:1:num_of_scales
    figure(i); t = tiledlayout(4,1);
    
    x = mra(i,1:K*l);
    x = reshape(x, K, l);

    % estimate higher order moments
    m   = sum(x,2) / l;
    var = sum((x-m).^2, 2) / l;
    std = sqrt(var);
    skw = sum((x-m).^3, 2) / l;
    krt = sum((x-m).^4, 2) / l;
    skw = skw ./ (std.^3);
    krt = krt ./ (std.^4);

    % subplot of reconstructed signal
    x1 = nexttile;
    plot(x1, time(1:K*l), sig1(1:K*l));
    xlabel('time in seconds'); 
    ylabel('Amplitude');
    
    % subplot of sliding variance
    x2 = nexttile;
    plot(x2, var);           
    ylabel('Variance');

    % subplot of sliding skewness
    x3 = nexttile;             
    plot(x3, skw);                 
    ylabel('Skewness');

    % subplot of sliding kurtosis
    x4 = nexttile;                     
    plot(x4, krt);           
    ylabel('Kurtosis');

    linkaxes([x1 x2 x3 x4], 'x');
end