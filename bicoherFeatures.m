% =======================================================
% Author Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------
%
% Function Description: 
%
% TODO: add text
% -------------------------------------------------------
%
% Arguments List: 
%
% TODO: add text
% -------------------------------------------------------
%
% Return List:
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

    % boolean masks for partitioning the bicoherence
    % matrices into delta waves, theta waves, alpha
    % waves and beta waves
    hex = freq + freq' + abs(freq + freq');
    deltaMask = (0  <= hex) & (hex < f1);
    thetaMask = (f1 <= hex) & (hex < f2);
    alphaMask = (f2 <= hex) & (hex < f3);
    betaMask  = (f3 <= hex) & (hex < f4);

    for i = 1:1:N
        % Extract the entire bicoherence
        % matrix for a 30sec EEG segment
        bic = cell2mat(X{i,1});

        % Partition the bicoherence matrix
        deltaBicoher = bic(deltaMask);
        thetaBicoher = bic(thetaMask);
        alphaBicoher = bic(alphaMask);
        betaBicoher  = bic(betaMask);
        
        % TODO: feature extraction from delta waves
        % TODO: feature extraction from theta waves
        % TODO: feature extraction from alpha waves
        % TODO: feature extraction from beta  waves
    end

    % Copy sleep stage Annotations;
    delta.Annotations = X.Annotations;
    theta.Annotations = X.Annotations;
    alpha.Annotations = X.Annotations;
    beta.Annotations  = X.Annotations;
end