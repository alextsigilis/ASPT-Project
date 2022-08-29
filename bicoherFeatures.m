% =======================================================
% Author Christodoulos Michaelides
% Date: August 25th, 2022
% -------------------------------------------------------
%
% Function Description: 
% This function can be used to extract features from the
% bicoherence estimations of EEG signals such as:
%   1) bicoherence entropy
%   2) bicoherence squared-entropy
%   4) average bicoherence
%   5) maximum bicoherence
% -------------------------------------------------------
%
% Arguments List: (X, freq)
% X: (table) a table with two columns. The first should
% contain the bicoherence estimations and the
% second should contain the sleep stage Annotations.
% You should use bicEEG to obtain this table.
% -------------------------------------------------------
%
% Return List: (delta, theta, alpha, beta)
%
% TODO: add text
% =======================================================

function [delta, theta, alpha, beta] = bicoherFeatures(X, f)
    % Number of EEG segments (30sec epochs)
    N = size(X,1);

    % Initialize empty tables to store
    % features from the bicoherence matrices.
    types = [
        "double",                   ...
        "double",                   ...
        "double",                   ...
        "double",                   ...
        "double",                   ...
        "string"];
    
    names = [                       ...
        "ent1",                     ...
        "ent2",                     ...
        "avg",                      ...
        "f1",                       ...
        "f2",                       ...
        "Annotations"];

    delta = table(                  ...
        'Size',          [N 6],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    theta = table(                  ...
        'Size',          [N 6],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    alpha = table(                  ...
        'Size',          [N 6],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    beta = table(                   ...
        'Size',          [N 6],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    % Binary masks for splitting the bicoherence matrix
    % into delta, theta, alpha and beta wave partitions
    f0 = 0.0; f1 = 4.0; f2 = 8.0; f3 = 16.0; f4 = 32.0;
    deltaHex = (f'<=f) & (f>=0) & (f0<=(f+f') & (f+f')<=f1);
    thetaHex = (f'<=f) & (f>=0) & (f1<=(f+f') & (f+f')<=f2);
    alphaHex = (f'<=f) & (f>=0) & (f2<=(f+f') & (f+f')<=f3);
    betaHex  = (f'<=f) & (f>=0) & (f3<=(f+f') & (f+f')<=f4);

    for i = 1:1:N
        % Extract the entire bicoherence
        % matrix for a 30sec EEG segment
        bic = cell2mat(X{i,1});

        % Partition the bicoherence matrix
        deltaBic = bic(deltaHex);
        thetaBic = bic(thetaHex);
        alphaBic = bic(alphaHex);
        betaBic  = bic(betaHex);
        
        % Extract features from delta waves
        p = deltaBic.^1; 
        q = deltaBic.^2; 

        delta{i,"ent1"} = -mean(p.*log2(p),'omitnan');
        delta{i,"ent2"} = -mean(q.*log2(q),'omitnan');
        delta{i,"f1"}   = nan;
        delta{i,"f2"}   = nan;
        delta{i,"avg"}  = mean(deltaBic);

        % Extract features from theta waves
        p = thetaBic.^1; 
        q = thetaBic.^2; 
        
        theta{i,"ent1"} = -mean(p.*log2(p),'omitnan');
        theta{i,"ent2"} = -mean(q.*log2(q),'omitnan');
        theta{i,"f1"}   = nan;
        theta{i,"f2"}   = nan;
        theta{i,"avg"}  = mean(thetaBic);

        % Extract features from alpha waves
        p = alphaBic.^1; 
        q = alphaBic.^2; 

        alpha{i,"ent1"} = -mean(p.*log2(p),'omitnan');
        alpha{i,"ent2"} = -mean(q.*log2(q),'omitnan');
        alpha{i,"f1"}   = nan;
        alpha{i,"f2"}   = nan;
        alpha{i,"avg"}  = mean(alphaBic);

        % Extract features from beta waves
        p = betaBic.^1; 
        q = betaBic.^2; 
        
        beta{i,"ent1"} = -mean(p.*log2(p),'omitnan');
        beta{i,"ent2"} = -mean(q.*log2(q),'omitnan');
        beta{i,"f1"}   = nan;
        beta{i,"f2"}   = nan;
        beta{i,"avg"}  = mean(betaBic);
    end

    % Copy sleep stage Annotations;
    delta.Annotations = X.Annotations;
    theta.Annotations = X.Annotations;
    alpha.Annotations = X.Annotations;
    beta.Annotations  = X.Annotations;
end