% =============================================================
% Author: Christodoulos Michaelides
% Date: August 22nd, 2022
% -------------------------------------------------------------
% 
% Demo Script 9
%
% 1) Choose a patient and an EEG recording from the dataset
% 2) Estimate the bicoherence of the EEG in every 30sec epoch
% 3) Plot the bicoherence estimations (pause for a few seconds
%    between each plot)
% =============================================================

% =============================================================
% 0) Reset your workspace variables
% =============================================================
clear all; close all; clc;

% =============================================================
% 1) Script parameters (Choose any value you want)
% =============================================================

patient_ID = 1;     % integer 1 ... 151
channel_ID = 1;     % integer 1 ... 4

% Parameters for bicoherence estimator
K  = 32;            % number of partitions per 30sec epoch
fs = 256;           % EEG sampling frequency in Hertz
fc = 32;            % cut-off frequency for bicoherence estimation

% Plot settings
dt = 1.0;           % Playback speed in seconds per frame     
levels = 8;         % number of levels for contour plots
cmap = "turbo";     % colormap for contour plots ("turbo" or "parula") 

% =============================================================
% 2) Bispectrum estimation
% =============================================================

fprintf("Loading EDF file ... "); tic;
Z = loadEDF(patient_ID);
fprintf("Done\n"); toc;

fprintf("Estimating bicoherence index ... "); tic;
[bic, freq] = bicEEG(Z,fs,fc,K,channel_ID);
fprintf("Done\n"); toc;

% =============================================================
% 3) Bicoherence Plots
% =============================================================

N = size(bic,1);
pause('on');

for i = 1:1:N
    b = cell2mat(bic{i,1});
    contourf(freq,freq,b,levels,'LineColor','none');
    txt = sprintf("frame %d out of %d: %s",i,N,bic.Annotations{i});
    xlabel('f_1'); ylabel('f_2'); title(txt);
    caxis([0.0 1.0]); colormap(cmap); colorbar;
    pause(dt);
end