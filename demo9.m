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
% 4) Save the contour plots as png images. (Use different 
%    folders based on the Sleep stage Annotations, EEG channel
%    and patient ID).
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
dt = 0.1;           % Playback speed in seconds per frame     
levels = 12;        % number of levels for contour plots
cmap = "turbo";     % colormap for contour plots ("turbo" or "parula") 

% save directory
dir = "C:\Users\USER\Desktop\bicoherence";

% =============================================================
% 2) Bicoherence estimation
% =============================================================

fprintf("Loading EDF file ... "); tic;
Z = loadEDF(patient_ID);
fprintf("Done\n"); toc;

fprintf("Estimating bicoherence index ... "); tic;
[bic, freq] = bicEEG(Z,K,fs,fc,channel_ID,"fancy");
fprintf("Done\n"); toc;

% =============================================================
% 3) Create empty folders to store the plots
% =============================================================

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\patient%d\\channel%d\\%s",                 ...
        dir, patient_ID, channel_ID, "Sleep stage W"));

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\patient%d\\channel%d\\%s",                 ...
        dir,patient_ID,channel_ID,"Sleep stage R"));

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\patient%d\\channel%d\\%s",                 ...
        dir,patient_ID,channel_ID,"Sleep stage N1"));

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\patient%d\\channel%d\\%s",                 ...
        dir,patient_ID,channel_ID,"Sleep stage N2"));

mkdir(                                                  ...
    sprintf(                                            ...
        "%s\\patient%d\\channel%d\\%s",                 ...
        dir,patient_ID,channel_ID,"Sleep stage N3"));

% =============================================================
% 4) Display contour plots and save the results
% =============================================================

N = size(bic,1);
pause('on');

for i = 1:1:N
    b = cell2mat(bic{i,1});
    
    contourf(freq,freq,b,levels,'LineColor','none');
    txt = sprintf("frame %d out of %d: %s",i,N,bic.Annotations{i});
    xlabel('f_1'); ylabel('f_2'); title(txt);
    caxis([0.0 1.0]); colormap(cmap); colorbar;
    
    path = sprintf("%s\\patient%d\\channel%d",dir,patient_ID,channel_ID);
    filename = sprintf("%s\\%s\\%d.png", path, bic.Annotations{i}, i);
    exportgraphics(gca, filename);
    
    pause(dt);
end