% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 19th, 2022
% -------------------------------------------------------------------
%
% Demo Script 7
%
% 1) Choose a patient from the dataset and an EEG channel 
%
% 2) Estimate the bispectrum of the EEG in a 30sec window
%
% 3) Display all the bispectra as seperate frames in a video 
% ===================================================================

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
% cmap: (string) colormap for contour plots (either "turbo" or "parula")
levels = 16;
dt = 1.0;
cmap = "turbo";

% ===================================================================
% 2) Estimate the Bispectra
% ===================================================================

fprintf("Loading EEG recordings ... "); tic;
X = loadEDF(patient_ID);
fprintf("Done\n"); toc;

fprintf("Estimating bispectra ... "); tic;
[bis, freq] = bisEEG(X,K,fs,fc,channel_ID,"fancy");
fprintf("Done\n"); toc;

% ===================================================================
% 3) Plot the bispectra. Pause for a few seconds between each plot
% ===================================================================

pause('on');

N = size(bis,1); 

for i = 1:1:N
    subplot(2,1,1);
    contourf(freq,freq,abs(cell2mat(bis{i,1})),levels,'LineColor','none');
    txt = sprintf("frame %d out of %d: %s",i,N,bis.Annotations{i});
    colormap(cmap); colorbar;
    xlabel("f_1"); ylabel("f_2"); title(txt);
    subplot(2,1,2);
    plot(cell2mat(X{i,channel_ID}));
    xlabel('samples'); ylabel('Amplitude');
    pause(dt);
end
