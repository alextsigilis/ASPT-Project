% =====================================================================
% Author: Christodoulos Michaelides
% Date: September 11th, 2022
% ---------------------------------------------------------------------
%
% Description: 
%
% 1) Select a number of patients and an EEG channel.
% 2) Estimate the bicoherence for every 30sec epoch.
% 3) Find the local maxima of the bicoherence matrices and  
%    detect QPC phenomena.
% 4) Plot histograms to find out whether or not there are any
%    notable changes in the distribution of the peaks depending
%    on the sleep stage

% Reset your Matlab workspace
clear all; close all; clc;

% Parameters for patient selection
first = 1;      % First patient ID
last  = 50;     % Last patient ID
channel = 1;    % Selected EEG channel

% Parameters for bicoherence estimation
K = 32;         % Number of segments
fs = 256;       % EEG sampling frequency
fc = 32;        % upper limit for bicoherence frequency axes

% Parameters for QPC detection
epsilon = 0.10; % Hard threshold for peak detection

% Create a table to store the local maxima of the bicoherence matrices
types = ["cell" "string"];
names = ["peaks" "Annotations"];
sz = [0 numel(names)];

QPC = table('Size',sz,'VariableTypes',types,'VariableNames',names);

for i = first:1:last

    % Make sure that the input files exist
    edf = sprintf("SN%03d.edf",i);
    mat = sprintf("%03d.mat",i);
    if (~isfile(edf)) && (~isfile(mat)) continue; end

    % Load the EEG recordings
    fprintf("Loading EDF files for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("OK\n");

    % Estimate the bicoherence of the selected EEG channel
    fprintf("Estimating bicoherence matrix ... ");
    [X,f] = bicEEG(Z,K,fs,fc,channel,"fast");
    fprintf("OK\n");

    % Extract features/peaks from the bicoherence
    fprintf("Locating bicoherence peaks ... ");
    QPC = [QPC; findQPC(X,f,epsilon)];
    fprintf("OK\n");

    fprintf("\n");
end


% Plot the distribution of the bicoherence peaks as a 2D-histogram
weightHist(QPC,f);