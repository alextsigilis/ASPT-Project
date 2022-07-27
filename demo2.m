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
% 0) Parameters for Wavelet Analysis
% ==================================================

input_file = "SN001.edf";        % name of input file for EEG recordings
channel = "EEGF4_M1";            % name of selected EEG/EOG/ECG channel
start = duration("00:00:00");    % starting point of EEG/EOG/ECG segment
range = duration("00:00:30");    % duration of EEG/EOG/ECG segment
wavelet = "db5";                 % wavelet type for MRA decomposition
fs = 256;                        % sampling frequency for EEG recording
levels = 7;                      % number of levels on MRA binary tree
pivot = 5;                       % Use every scale after pivot for 
                                 % reconstruction (including pivot)

% =========================================================================
% Do NOT change anything below that point,
% unless you know what you are doing.

% Construct the time axis of the EEG segment
t0 = seconds(start);
t1 = seconds(start+range);
time = linspace(t0,t1,fs*(t1-t0));

num_of_scales = levels + 1;        % number of frequency scales
scales = (fs/2) * 2.^-[0:levels];  % frequency boundaries for every scale
scales = [scales 0];

% ==================================================
% 1) Read EEG recordings from the EDF file
% ==================================================

 % Read the EDF file
X = edfread(input_file);

% extract a small segment from the selected EEG channel
interval = timerange(start,start+range);
sig = X{interval,channel};      

% Convert cell-array into 1D-array
sig = cell2mat(sig);          

% ==================================================
% 2) Perform MRA Decomposition and Reconstruction
% ==================================================

% Logical array for selecting reconstruction elements
levelForReconstruction = (1:1:levels) >= pivot;

% Perform the decomposition using modwt
wt = modwt(sig,wavelet,levels);

% Construct MRA matrix using modwtmra
mra = modwtmra(wt,wavelet);

% Sum down the rows of the selected multiresolution signals
sig1 = sum(mra(levelForReconstruction,:),1);

% ==========================================================
% 3) Plot the results
% ==========================================================

% Original (sig) vs Reconstructed (sig1)
figure(1); 
plt = plot(time,sig, time,sig1);
plt(1).LineWidth = 0.5;
plt(2).LineWidth = 2;
xlabel('time in seconds');
ylabel('Amplitude in microvolts');
title('Original vs Reconstructed Signal');
legend("original", "reconstructed");

% Plot every frequency scale (both selected and discarded ones)
figure(2);

for i = 1:1:num_of_scales
    subplot(num_of_scales,1,i);
    plot(time, mra(i,:));
    xlabel('time');
    ylabel('Amplitude');
    title('scale ' + string(i));
end

% Plot a histogram for every frequency scale
figure(3);

for i = 1:1:num_of_scales
    subplot(num_of_scales,1,i);
    hist(mra(i,:),30);
    title('scale ' + string(i));
end

% =========================================================
% 4) 2nd, 3rd and 4th order statistics on all frequency
% scales
% =========================================================

var = std(mra,0,2);       % second order statistics on all scales
skw = skewness(mra,0,2);  % third order statistics on all scales
krt = kurtosis(mra,0,2);  % fourth order statistics on all scales

for i = 1:1:num_of_scales
    % frequency boundaries of current scale
    f1 = scales(i+1); f2 = scales(i);
	
    fprintf("Frequency scale: %dHz-%dHz\n",f1,f2);
    fprintf("variance = %f\n", var(i));
	fprintf("skewness = %f\n", skw(i));
    fprintf("kurtosis = %f\n", krt(i));
    fprintf("==========================\n");
    fprintf("\n");
end