clear all; close all; clc;

% ------------------------ Script Parameters ------------------------

% Hyperparameters for estimating the bicoherence
K  = 24;                        % Number of segments
fs = 256;                       % Sampling frequency
fc = 32;                        % upper bound on frequency axis            

% Select patients and an EEG channel from the dataset
channel = 4;                    % EEG channel
start = 1;                      % first patient
stop = 50;                      % last patient

% Plot settings
nbins = [8 30 30 30 30];    % Number of histogram bins for each feature

% ---------------- Do not change anything below ---------------------

% Initialize empty tables to store bicoherence features
names = [           ...
    "numOfPeaks",   ...
    "avgBic",       ...
    "avgDist",      ...
    "maxDist",      ...
    "center",       ...
    "Annotations"];

types = [repmat("double",1,numel(names)-1) "string"];
sz = [0 numel(names)];

Y = table(                      ...
    'Size',          sz,        ...
    'VariableTypes', types,     ...
    'VariableNames', names);


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

    % Locate the peaks of the bicoherence matrix
    fprintf("Locating bicoherence peaks ... ");
    QPC = findQPC(X);
    fprintf("OK\n");

    % Extract features from every region of the 
    % bicoherence matrices
    fprintf("Extracting QPC features ... ");
    Y = [Y; QPCfeatures(QPC,f)];
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
    [y1, x1] = hist(z1,nbins(k)); y1 = (y1 * nbins(k)) ./ (numel(z1) * range(z1));
    [y2, x2] = hist(z2,nbins(k)); y2 = (y2 * nbins(k)) ./ (numel(z2) * range(z2));
    [y3, x3] = hist(z3,nbins(k)); y3 = (y3 * nbins(k)) ./ (numel(z3) * range(z3));
    [y4, x4] = hist(z4,nbins(k)); y4 = (y4 * nbins(k)) ./ (numel(z4) * range(z4));
    [y5, x5] = hist(z5,nbins(k)); y5 = (y5 * nbins(k)) ./ (numel(z5) * range(z5));
 
    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(Y.Properties.VariableNames{k}); 
    ylabel('Probability Density Function');
    legend("W", "R", "N1", "N2", "N3");
end