% =======================================================
% Author Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------
%
% Function Description: 
% This function can be used to extract features from the
% bispectrum and bicoherence estimations of EEG signals
% such as:
% => Entropy of bicoherence
% => Entropy of squared bicoherence/bispectrum
% => Average of log-bispectrum/log-bicoherence
% => Average of main diagonal 
% => 1st order spectral moment of main diagonal
% -------------------------------------------------------
%
% Arguments List: (X, freq)
% X: (table) a table with two columns. The first should
% contain the bispectrum/bicoherence estimations and the
% second should contain the sleep stage Annotations.
% You should use bisEEG or bicEEG to obtain this table.
% 
% freq (1D vector) the frequency axis of the
% bicoherence/bispectrum matrices. You should use 
% bisEEG or bicEEG to obtain this vector.
% -------------------------------------------------------
%
% Return List: (delta, theta, alpha, beta)
%
% TODO: add text
% =======================================================

function [delta, theta, alpha, beta] = bicoherFeatures(X,freq)
    % Number of EEG segments (30sec epochs)
    N = size(X,1);

    % f1: (float) maximum frequency of delta waves
    % f2: (float) maximum frequency of theta waves
    % f3: (float) maximum frequency of alpha waves
    % f4: (float) maximum frequency of beta  waves
    f1 = 4.0; f2 = 8.0; f3 = 16.0; f4 = 32.0;

    % Initialize empty tables to store features
    % from the bicoherence / bispectrum matrices.
    % delta: features from delta waves
    % theta: features from theta waves
    % alpha: features from alpha waves
    % beta:  features from beta waves

    types = [           ...
        "double",       ...
        "double",       ...
        "double",       ...
        "double",       ...
        "double",       ...
        "string"];

    names = [           ...
        "ent",          ...     % entropy
        "entsqr",       ...     % squared entropy
        "logsum",       ...     % sum of logarithmic bicoherence
        "logsumdiag",   ...     % diagonal sum of logarithmic bicoherence
        "logmoment",    ...     % 1st order spectral moment of main diagonal
        "Annotations"];

    delta = table(                  ...
        'Size',          [N 5],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    theta = table(                  ...
        'Size',         [N 5],      ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    alpha = table(                  ...
        'Size',          [N 5],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    beta  = table(                  ...
        'Size',          [N 5],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    % boolean masks for partitioning the bicoherence
    % matrices into delta waves, theta waves, alpha
    % waves and beta waves
    hex = abs(freq) + abs(freq') + abs(freq + freq');
    deltaMask =    (0 <= hex) & (hex < 2*f1);
    thetaMask = (2*f1 <= hex) & (hex < 2*f2);
    alphaMask = (2*f2 <= hex) & (hex < 2*f3);
    betaMask  = (2*f3 <= hex) & (hex < 2*f4);

    for i = 1:1:N
        % Extract the entire bicoherence / bispectrum
        % matrix for a 30sec EEG segment
        bic = cell2mat(X{i,1});

        % Partition the bicoherence / bispectrum matrix
        deltaBicoher = bic(deltaMask);
        thetaBicoher = bic(thetaMask);
        alphaBicoher = bic(alphaMask);
        betaBicoher  = bic(betaMask);
        
        % TODO: feature extraction from delta waves
        delta{i,"ent"}        = nan;
        delta{i,"entsqr"}     = nan;
        delta{i,"logsum"}     = nan;
        delta{i,"logsumdiag"} = nan;
        delta{i,"logmoment"}  = nan;

        % TODO: feature extraction from theta waves
        theta{i,"ent"}        = nan;
        theta{i,"entsqr"}     = nan;
        theta{i,"logsum"}     = nan;
        theta{i,"logsumdiag"} = nan;
        theta{i,"logmoment"}  = nan;

        % TODO: feature extraction from alpha waves
        alpha{i,"ent"}        = nan;
        alpha{i,"entsqr"}     = nan;
        alpha{i,"logsum"}     = nan;
        alpha{i,"logsumdiag"} = nan;
        alpha{i,"logmoment"}  = nan;

        % TODO: feature extraction from beta  waves
        beta{i,"ent"}        = nan;
        beta{i,"entsqr"}     = nan;
        beta{i,"logsum"}     = nan;
        beta{i,"logsumdiag"} = nan;
        beta{i,"logmoment"}  = nan;
    end

    % Copy sleep stage Annotations;
    delta.Annotations = X.Annotations;
    theta.Annotations = X.Annotations;
    alpha.Annotations = X.Annotations;
    beta.Annotations  = X.Annotations;
end