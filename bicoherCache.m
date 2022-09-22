% ========================================================================
% Author: Christodoulos Michaelides
% Date: September 21st, 2022
% ------------------------------------------------------------------------
%
% Description: 
% Estimating the bicoherence matrices for every patient, every EEG channel
% and every 30sec epoch is an extremely time-consuming workload even when
% using an efficient algorithm for estimating those matrices. This script
% estimates those matrices for a user-defined number of patients and saves 
% those results in a mat-file for future use (eg feature extraction).
% ------------------------------------------------------------------------
%
% Script parameters:
% 
% start: (integer) first patient selected from the dataset
% stop:  (integer) last patient selected from the dataset
% cacheFile: (string) name of the final mat file.
% K:      (integer) number of segments when estimating the bicoherence
% fs:     (integer) EEG-sampling frequency in Hertz
% fc:     (integer) upper-bound on the frequency axes in Hertz
% method: (string) estimation technique for bicoherence matrices
%         ==> "fancy" estimates the bicoherence on the entire 2D plane
%         ==> "fast" estimates the bicoherence on its primary region only
% ------------------------------------------------------------------------
% 
% Output File:
% cacheFile: (string) the name of the output file where the following
% workspace variables are stored:
%   => f:           (1D array of floats) the frequeny axis in Hertz
%   => K:           See previous paragraph
%   => fs:          See previous paragraph
%   => fc:          See previous paragraph
%   => method:      See previous paragraph
%   => bicTable:    a table with 4+1 columns. The first 4 columns contain
%                   the bicoherence matrices for every EEG channel and 
%                   every 30sec epoch (estimated with the user-defined
%                   parameters). The 5th column contains the Sleep stage
%                   Annotations.
% ========================================================================

% Reset your Matlab workspace
clear all; close all; clc;

% Script parameters (See documentation)
start = 1;
stop  = 50;

cacheFile = "bicoherCache.mat";

K  = 32;
fs = 256;
fc = 32;
method = "fast";

% --------------- Do not change anything below that point ---------------

% Delete previous cache files
if isfile(cacheFile)
    fprintf("Removing previous cache files ... ");
    delete(cacheFile);
    fprintf("Done\n\n");
end

% Initialise an empty table to store the bicoherence matrices and the 
% sleep stage Annotations
sz    = [0 5];
types = ["cell" "cell" "cell" "cell" "string"];
names = ["EEGF4_M1" "EEGC4_M1" "EEGO2_M1" "EEGC3_M2" "Annotations"];

bicTable = table(           ...
    'Size',sz,              ...
    'VariableTypes',types,  ...
    'VariableNames',names);

% Iterate through all the selected patients
for i = start:1:stop
    % Make sure that the input file exists
    edfFile = sprintf("SN%03d.edf",i);
    matFile = sprintf("%03d.mat",i);

    if (~isfile(edfFile)) && (~isfile(matFile))
        continue;
    end

    % Progress status
    tic; 
    fprintf("Patient: %d\n", i);

    % Load the EEG recordings
    fprintf("Loading EEG recordings\n");
    Z = loadEDF(i); N = size(Z,1);
    
    % bicoherence on the first EEG channel
    fprintf("Estimating bicoherence in EEG channel 1\n");
    [b1, ~] = bicEEG(Z,K,fs,fc,1,method);

    % bicoherence on the second EEG channel
    fprintf("Estimating bicoherence in EEG channel 2\n");
    [b2, ~] = bicEEG(Z,K,fs,fc,2,method);

    % bicoherence on the third EEG channel
    fprintf("Estimating bicoherence in EEG channel 3\n");
    [b3, ~] = bicEEG(Z,K,fs,fc,3,method);

    % bicoherence on the fourth EEG channel
    fprintf("Estimating bicoherence in EEG channel 4\n");
    [b4, f] = bicEEG(Z,K,fs,fc,4,method);

    % Save the bicoherence matrices in the table
    fprintf("Saving results in memory\n");
    temp = table(               ...
        'Size',[N 5],           ...
        'VariableTypes',types,  ...
        'VariableNames',names);
    
    temp{:,1} = b1{:,1};
    temp{:,2} = b2{:,1};
    temp{:,3} = b3{:,1};
    temp{:,4} = b4{:,1};
    temp{:,5} = b4{:,2};
    
    bicTable = [bicTable; temp];

    toc;

    % Progress status
    fprintf("\n");
end

% Copy the bicoherence matrices in a cache file (mat file)
fprintf("Memory dump. Please wait ... ");
save(cacheFile,"bicTable", "f", "K", "fs", "fc", "method");
fprintf("Done\n\n");

% Progress status
fprintf("Cache file directory: %s\n", pwd);
fprintf("Cache file name: %s\n\n", cacheFile);

% Reset your workspace
clear all; close all;