
% ==========================================================
% Reset your workspace variables
% ==========================================================

clear all; close all; clc;

% ==========================================================
% Script parameters (Choose any value you want)
% ==========================================================

id1 = 1;            % first patient
id2 = 25;           % last patient

nbins = 100;        % number of bins for histograms

% ==========================================================
% Extract features from EOG recordings
% ==========================================================

features = table(                                   ...
    'Size',             [0 2],                      ...
    'VariableTypes',    ["double" "string"],        ...
    'VariableNames',    ["xcorr" "Annotations"]);

for i = id1:1:id2
    if ~isfile(sprintf("SN%03d.edf",i)) continue; end

    fprintf("Loading EOG recordings for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("Done\n");

    fprintf("Extracting features from EOG recordings ... ");
    features = [features; featuresEOG(Z)];
    fprintf("Done\n");
end

% ==========================================================
% Plot histograms/scatterplots of features
% ==========================================================

% mask1: (1D array) boolean array for REM epochs
% mask2: (1D array) boolean array for NREM epochs
mask1 = features{:,"Annotations"} == "Sleep stage R";
mask2 = features{:,"Annotations"} ~= "Sleep stage R";

% Histogram of cross-correlation for REM sleep 
[counts1, centers1] = hist(features{mask1, "xcorr"}, nbins);
counts1 = counts1 / sum(counts1);

% Histogram of cross-correlation for NREM sleep 
[counts2, centers2] = hist(features{mask2, "xcorr"}, nbins);
counts2 = counts2 / sum(counts2);

% Display the histograms
figure(1); hold on; grid on;
plot(centers1,counts1,'r',centers2,counts2,'b');
xlabel("EOG cross-correlation coefficients");
ylabel("Relative frequency");
title("Histogram of cross-correlation");
legend("REM","Non-REM");