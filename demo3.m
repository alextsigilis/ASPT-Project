% ======================================================
% WARNING: Don't even bother running this script.
% There are many bugs that need to be fixed.
% ======================================================

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
% to reconstruct an approximation of the orignal
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

clear all; clc;     % reset your workspace and clear terminal

% ======================================================
% 1) Parameters of MRA Decomposition and Reconstruction
% (Choose any value that you want to experiment with)
% ======================================================

input_file = "SN001.edf";           % name of input file
channel = "EEGF4_M1";               % name of selected channel
win_length = duration("00:00:03");  % window length when estimating moments
wavelet = "db5";                    % wavelet for MRA 
levels = 7;                         % number of levels on MRA binary tree 
fs = 256;                           % sampling frequency of the recording
pivot = 5;                          % Use every frequency scale below pivot  
                                    % (including the pivot scale) for 
                                    % signal reconstruction

% ------------------------------------------------------
% Do not change anything below that point, 
% unless you know exactly what you are doing.
% ------------------------------------------------------

% number of frequency scales
num_of_scales = levels + 1;                

% frequency boundaries for every scale
scales = (fs/2) * 2 .^ 1:1:num_of_scales;
scales = [scales 0];

% time axis for EEG/EOG/ECG signal
time = linspace(0, ...);             % TODO: fix that shit

% ======================================================
% 2) Read data from the EDF file
% ======================================================

X = edfread("SN001.edf");        % Read the entire file as a timetable
sig = X{:,channel};              % choose the appropriate channel
sig = cell2mat(sig);             % convert chosen channel to 1D-array

clear X;                         % delete the timetable to 
                                 % save some space in RAM

% ======================================================
% 3) Perform MRA decomposition and reconstruction
% ======================================================

% Logical array for selecting reconstruction scales
levelForReconstruction = (1:1:num_of_scales) >= pivot;

% Perform the decomposition using modwt
wt = modwt(sig,wavelet,levels);

% Construct MRA matrix using modwtmra
mra = modwtmra(wt,wavelet);

% Sum down the rows of the selected multiresolution signals
sig1 = sum(mra(levelForReconstruction,:),1);

% Delete unused variables to save some space on RAM
clear wt;  

% ======================================================
% 4) Estimate variance, skewness and kurtosis in a 
% sliding window 
% ======================================================

% Convert the reconstructed signal 
% into a 2D array.
% ---------------------------------------
% Every row of this 2D array contains
% a small window  of the reconstructed
% signal where variance, skewness and
% kurtosis are locally estimated.
% ---------------------------------------
% The duration of the window is adjusted
% in section 1) of the script 

n = fs * win_length;            % size of sliding window in samples

x = reshape(sig1, ...);         % TODO: fix that shit

m   = sum(x, 1) / n;            % sliding average
var = sum((x-m).^2, 1) / n;     % sliding variance (biased estimate) 
std = sqrt(var);                % sliding standard deviation

skw = sum((x-m).^3, 1) / n;     % sliding 3rd moment
krt = sum((x-m).^4, 1) / n;     % sliding 4th moment

skw = skw ./ (std.^3);          % sliding skewness (biased estimate)
krt = krt ./ (std.^4);          % sliding kurtosis (biased estimate)

% delete unused variables
clear x, avg;                 

% ======================================================
% 5) Plot the results
% ======================================================

% Reconstructed vs Original Signal
figure(1); 

plt = plot(time,sig, time, sig2);
plt(1).LineWidth = 0.5; 
plt(2).LineWidth = 2.0;

xlabel("time in seconds");
ylabel("Amplitude in microVolts");
title("Original vs Reconstruction");

% ------------------------------------------------------

% Plots of selected frequency scales
% and sliding window statistics
for i = pivot:1:num_of_scales
    figure(i);

    % upsample std(i,:), var(i,:) skw(i,:) 
    % and krt(i,:) in order to match the 
    % length of the reconstructed signal
    std = interp1(std(i,:), time);
    var = interp1(var(i,:), time);
    skw = interp1(skw(i,:), time);
    krt = interp1(krt(i,:), time);

    % subplot of reconstructed signal
    subplot(5,1,1);
    plot(time, sig1);
    xlabel('time in seconds'); 
    ylabel('Amplitude');
    
    % subplot of frequency scale
    subplot(5,1,2);        
    plot(time, mra(i,:));
    xlabel('time in seconds'); 
    ylabel('Amplitude');
    
    % subplot of sliding variance
    subplot(5,1,3);
    plot(time, var);           
    xlabel('time in seconds'); 
    ylabel('Variance');

    % subplot of sliding skewness
    subplot(5,1,4);              
    plot(time, skw);                
    xlabel('time in seconds'); 
    ylabel('Skewness');

    % subplot of sliding kurtosis
    subplot(5,1,5);                      
    plot(time, krt);           
    xlabel('time in seconds'); 
    ylabel('Kurtosis');
end