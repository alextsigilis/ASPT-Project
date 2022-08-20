% ===================================================================
% 0) Clear your workspace variables
% ===================================================================

clear all; close all; clc;

% ===================================================================
% 1) Script Parameters (choose any value you want to experiment with)
% ===================================================================

% patient_ID: (integer) 1 ... 151
% channel_ID: (integer) 1 ... 4
patient_ID = 1;
channel_ID = 1;

% K:  (integer) number of partitions per 30sec EEG recording
% fs: (float) sampling frequency for EEG recordings
% fc: (float) truncation frequency for bispectrum estimations
K  = 4;
fs = 256;
fc = 4;

% levels: (integer) number of levels for contour plots
% dt: (integer) playback speed in seconds per frame
levels = 8;
dt = 1.0;

% ===================================================================
% 2) Estimate the Bispectra
% ===================================================================

fprintf("Loading EEG recordings ... "); tic;
X = loadEDF(patient_ID);
fprintf("Done\n"); toc;

fprintf("Estimating bispectra ... "); tic;
[bis, freq] = bisEEG(X,K,fs,fc,channel_ID);
fprintf("Done\n"); toc;

clear X;

% ===================================================================
% 3) Plot the bispectra. Pause for a few seconds between each plot
% ===================================================================

pause('on');

N = size(bis,1); 

for i = 1:1:N
    contourf(freq, freq, abs(cell2mat(bis{i,1})), levels);
    txt = sprintf("frame %d out of %d: %s", i, N, bis.Annotations{i});
    xlabel("f_1"); ylabel("f_2"); title(txt);
    pause(dt);
end