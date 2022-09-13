% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------------------
%
% Description:  
%
% This script estimates the EEG-bicoherence matrices on a large 
% number of patients. It then extracts features from the bicoherence 
% estimates and plots histograms of those features. Different
% histograms are generated for every sleep stage. We can identify
% useful features by visually inspecting those histograms.
% ===================================================================

% Clear all workspace variables,
% Close any open windows and
% clear the command line.
clear all; close all; clc;

% ------------------------ Script Parameters ------------------------

% Hyperparameters for estimating the bicoherence
K  = 15; K = 32;                % Number of segments
fs = 256;                       % Sampling frequency
fc = 32;                        % upper bound on frequency axis

% Select patients and an EEG channel from the dataset
channel = 4;                    % EEG channel
start = 1;                      % first patient
stop = 50;                      % last patient

% Plot settings
nbins = 100;                    % Number of histogram bins
norm = "pdf";                   % Normalization type for histograms
                                % cdf => cumulative distribution 
                                % pdf => probability density function

% ---------------- Do not change anything below ---------------------

% Initialize empty tables to store bicoherence features
types = ["double" "double" "double" "string"];
names = ["ent1" "ent2" "avg" "Annotations"];
sz = [0 numel(names)];

Y = table('Size',sz,'VariableTypes',types,'VariableNames',names);

% Extract bicoherence features from every patient
for i = start:stop
    % Make sure the input file exists
    edf = sprintf("SN%03d.edf",i);
    mat = sprintf("%03d.mat",i);
    if (~isfile(edf)) && (~isfile(mat)) continue; end

    % Load the EEG recordings 
    fprintf("Loading EDF files for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("OK\n");

    % Choose an EEG channel and estimate the bicoherence
    fprintf("Estimating bicoherence matrix ... ");
    [X,f] = bicEEG(Z,K,fs,fc,channel,"fast");
    fprintf("OK\n");

    % Extract features from every region of the 
    % bicoherence matrices
    fprintf("Extracting bicoherence features ... ");
    Y = [Y; bicoherFeatures(X,f)];
    % [D, U, A, B] = bicoherFeatures(X,f);
    fprintf("OK\n");

    fprintf("\n");
end

% ----------------- Histograms of bicoherence features -----------------

% incremental counter to distinguish between different plots and figures
idx = 1;

% binary arrays to split features based on sleep stage annotations
W  = Y.Annotations == "Sleep stage W";
R  = Y.Annotations == "Sleep stage R";
N1 = Y.Annotations == "Sleep stage N1";
N2 = Y.Annotations == "Sleep stage N2";
N3 = Y.Annotations == "Sleep stage N3";

K = size(Y,2) - 1;

% Plot histograms for every feature
for k = 1:K
    z1 = Y{W,k};        % Features of sleep stage W
    z2 = Y{R,k};        % Features of sleep stage R
    z3 = Y{N1,k};       % Features of sleep stage N1
    z4 = Y{N2,k};       % Features of sleep stage N2
    z5 = Y{N3,k};       % Features of sleep stage N3

    % Probability Density Functions for different sleep stages
    [y1, x1] = hist(z1,nbins); y1 = (y1 * nbins) ./ (numel(z1) * range(z1));
    [y2, x2] = hist(z2,nbins); y2 = (y2 * nbins) ./ (numel(z2) * range(z2));
    [y3, x3] = hist(z3,nbins); y3 = (y3 * nbins) ./ (numel(z3) * range(z3));
    [y4, x4] = hist(z4,nbins); y4 = (y4 * nbins) ./ (numel(z4) * range(z4));
    [y5, x5] = hist(z5,nbins); y5 = (y5 * nbins) ./ (numel(z5) * range(z5));
 
    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(Y.Properties.VariableNames{k}); 
    ylabel('Probability Density Function');
    legend("W", "R", "N1", "N2", "N3");
end