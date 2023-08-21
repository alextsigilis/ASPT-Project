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

cacheFile = "bicoherCache.mat";

K  = 32;
fs = 256;
fc = 32;
method = "fast";

% --------------- Do not change anything below that point ---------------

% Disable all warning
warning('off','all');

% Select edf files
[files, path] = uigetfile("*.edf", 'MultiSelect', 'on');

% Delete previous cache files
if isfile(cacheFile)
    fprintf("Removing previous cache files ... ");
    delete(cacheFile);
    fprintf("Done\n\n");
end

% Initialise an empty table to store the bicoherence matrices and the 
% sleep stage Annotations
nchan = 8;
sz    = [0 9];
types = ["cell" "cell" "cell" "cell" "cell" "cell" "cell" "cell" "string"];
names = ["EEGF4_M1", "EEGC4_M1", "EEGO2_M1", "EEGC3_M2", "EMGChin", "EOGE1_M2", "EOGE2_M2", "ECG", "Annotations"];


bicTable = table(           ...
    'Size',sz,              ...
    'VariableTypes',types,  ...
    'VariableNames',names);

% Iterate through all the selected patients
for i = 1:numel(files)
    
    file = files{i};
    pid = str2num(file(3:5));

    % Make sure that the input file exists
    edfFile = fullfile(path, sprintf("SN%03d.edf",pid));
    matFile = sprintf("%03d.mat",pid);

    if (~isfile(edfFile)) && (~isfile(matFile))
        continue;
    end

    % Progress status
    tic; 
    fprintf("Patient: %d\n", pid);
    
    % Load the EEG recordings
    fprintf("Loading EEG recordings\n");
    Z = loadEDF(pid, path); N = size(Z,1);
    
    b = {};
    
    for j = 1:nchan
        % bicoherence on the j-th channel
        fprintf("Estimating bicoherence in channel %d\n", j);
        [b{j}, ~] = bicEEG(Z,K,fs,fc,1,method);
    end

    % Save the bicoherence matrices in the table
    fprintf("Saving results in memory\n");
    temp = table(               ...
        'Size',[N sz(2)],           ...
        'VariableTypes',types,  ...
        'VariableNames',names);
    
    for j = 1:nchan
        bj = b{j};
        temp{:,j} = bj{:,1};
    end
    temp{:,9} = bj{:,2}; % !!!! The table has 9 collumns, the last on is the annotations.
    
    bicTable = [bicTable; temp];

    toc;

    % Progress status
    fprintf("\n");
end

% Copy the bicoherence matrices in a cache file (mat file)
fprintf("Memory dump. Please wait ... ");
save(cacheFile,"bicTable", "K", "fs", "fc", "method");
fprintf("Done\n\n");

% Progress status
fprintf("Cache file directory: %s\n", pwd);
fprintf("Cache file name: %s\n\n", cacheFile);

% Reset your workspace
clear all; close all;