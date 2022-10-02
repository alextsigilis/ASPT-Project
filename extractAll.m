% ===================================================================
% Authors:    Chrysa Doulou
%             Christodoulos Michaelides
%             Demetrios Orphanos
%             Stergios Gregoriou
%
% Date:       October 2nd, 2022
% -------------------------------------------------------------------
%
% Script Description:
% This script extracts a wide variety of useful features from the
% polysomnographic recordings of the dataset. Those features include:
%   1) Bicoherence Features
%   2) Cepstrum Features
%   3) DWT features
%
% Those features are stored in .mat files and can later be used to 
% to train classifiers for automatic sleep stage scoring.
% -------------------------------------------------------------------

% ===================================================================
% Clear all previous workspace variables
% ===================================================================

clear all;      % Delete workspace variables 
close all;      % Close all open windows
clc;            % Clear the terminal

% ===================================================================
% Script Parameters
% ===================================================================

first = 1;          % First patient selected from the dataset
last  = 20;         % Last patient selected from the dataset

fs = 256;           % Sampling frequency of PSG recordings (in Hertz)
dt = 30;            % Epoch duration of PSG recordings (in seconds)

K  = 24;            % Number of partitions for bicoherence matrices
fc = 32;            % Maximum frequency for bicoherence matrices

% Save folder for feature files. (-mat files)
path = sprintf("C:\\Users\\USER\\Desktop\\features");

% ===================================================================
% Feature Extraction
% ===================================================================

for n = first:1:last

    % Make sure that the EDF/mat file exists
    edfFile = sprintf("SN%03d.edf",n);
    matFile = sprintf("%03d.mat",n);

    if (~isfile(edfFile)) && (isfile(matFile))
        continue;
    end

    % Progress status
    fprintf("Patient %d:\n",n);

    % Load PSG recordings
    fprintf("Loading PSG recordings ... ");
    Z = loadEDF(n);
    fprintf("Done\n");

    % Preprocessing
    fprintf("Prefiltering ... ");
    Z = prefilter(Z,fs);
    fprintf("Done\n");

    % ---------------- DWT Features ----------------
    
    fprintf("Performing MRA decomposition ... ");
    coeff1 = mraEEG(Z,"EEGF4_M1");
    coeff2 = mraEEG(Z,"EEGC4_M1");
    coeff3 = mraEEG(Z,"EEGO2_M1");
    coeff4 = mraEEG(Z,"EEGC3_M2");
    fprintf("Done\n");

    fprintf("Extracting features from DWT coefficients ... ");
    [D1, T1, A1, B1] = statEEG(coeff1);
    [D2, T2, A2, B2] = statEEG(coeff2);
    [D3, T3, A3, B3] = statEEG(coeff3);
    [D4, T4, A4, B4] = statEEG(coeff4);
    fprintf("Done\n");

    % TODO: store features in tables

    % ------------ Bicoherence Features ------------

    fprintf("Estimating bicoherence matrices ... ");
    [b1, f1] = bicEEG(Z,K,fs,fc,"EEGF4_M1","fast");
    [b2, f2] = bicEEG(Z,K,fs,fc,"EEGC4_M1","fast");
    [b3, f3] = bicEEG(Z,K,fs,fc,"EEGO2_M1","fast");
    [b4, f4] = bicEEG(Z,K,fs,fc,"EEGC3_M2","fast");
    fprintf("Done\n");

    fprintf("Extracting Features from bicoherence matrices ... ");
    Y1 = bicoherFeatures(b1,f1);
    Y2 = bicoherFeatures(b2,f2);
    Y3 = bicoherFeatures(b3,f3);
    Y4 = bicoherFeatures(b4,f4);
    fprintf("Done\n");

    % TODO: store features in tables

    % ---------------- QPC features ----------------
    fprintf("Locating QPC frequency-pairs ... ");
    qpc1 = findQPC(b1);
    qpc2 = findQPC(b2);
    qpc3 = findQPC(b3);
    qpc4 = findQPC(b4);
    fprintf("Done\n");

    fprintf("Extracting QPC features ... ");
    Y1 = QPCfeatures(qpc1,f1);
    Y2 = QPCfeatures(qpc2,f2);
    Y3 = QPCfeatures(qpc3,f3);
    Y4 = QPCfeatures(qpc4,f4);
    fprintf("Done\n");

    % TODO: store features in tables

    % ------------- Cepstrum Features --------------

    % ---------------- EOG features ----------------

    % ---------------- EMG features ----------------
end