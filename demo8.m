% ==========================================================
% Reset your workspace variables
% ==========================================================

clear all; close all; clc;

% ==========================================================
% Script parameters (Choose any value you want)
% ==========================================================

id1 = 1;            % first patient
id2 = 50;           % last patient
nbins = 100;        % number of bins for histograms

% ==========================================================
% Extract features from EOG recordings
% ==========================================================

sz = [0 6];
types = ["double" "double" "double" "double" "double" "string"];
names = ["xcorr" "stdDev1" "stdDev2" "AUC1" "AUC2" "Annotations"];
features = table('Size',sz,'VariableTypes',types,'VariableNames',names);

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
% mask3: (1D array) boolean array for wakefullness (Sleep stage W) 
mask1 = features{:,"Annotations"} == "Sleep stage R";
mask2 = features{:,"Annotations"} ~= "Sleep stage R";
mask3 = features{:,"Annotations"} ~= "Sleep stage W";

% ----------------------------------------------------------------

% Histogram of cross-correlation for REM sleep
[y1, x1] = hist(features{mask1 & mask3, "xcorr"}, nbins);
y1 = y1 / sum(y1);

% Histogram of cross-correlation for NREM sleep 
[y2, x2] = hist(features{mask2 & mask3, "xcorr"}, nbins);
y2 = y2 / sum(y2);

% Display the histograms
figure(1); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b');
xlabel("EOG cross-correlation coefficients");
ylabel("Relative frequency");
title("Histogram of EOG cross-correlation");
legend("REM","Non-REM");

% ----------------------------------------------------------------

% Histogram of AUC (EOG E1-M2) for REM sleep
[y1, x1] = hist(features{mask1 & mask3, "AUC1"}, nbins);
y1 = y1 / sum(y1);

% Histogram of AUC (EOG E1-M2) for NREM sleep
[y2, x2] = hist(features{mask2 & mask3, "AUC1"}, nbins);
y2 = y2 / sum(y2);

% Display the Histograms
figure(2); hold on; grid on; 
plot(x1,y1,'r',x2,y2,'b');
xlabel("AUC values");
ylabel("Relative frequency");
title("Histogram of AUC for EOG E1-M2")
legend("REM", "Non-REM");

% ----------------------------------------------------------------

% Histogram of AUC (EOG E2-M2) for REM sleep
[y1, x1] = hist(features{mask1 & mask3, "AUC2"}, nbins);
y1 = y1 / sum(y1);

% Histogram of AUC (EOG E2-M2) for NREM sleep
[y2, x2] = hist(features{mask2 & mask3, "AUC2"}, nbins);
y2 = y2 / sum(y2);

% Display the Histograms
figure(3); hold on; grid on; 
plot(x1,y1,'r',x2,y2,'b');
xlabel("AUC values");
ylabel("Relative frequency");
title("Histogram of AUC for EOG E2-M2")
legend("REM", "Non-REM");

% ----------------------------------------------------------------

% Histogram of standard deviation (EOG E1-M2) for REM sleep
[y1, x1] = hist(features{mask1 & mask3, "stdDev1"}, nbins);
y1 = y1 / sum(y1);

% Histogram of standard deviation (EOG E1-M2) for NREM sleep
[y2, x2] = hist(features{mask2 & mask3, "stdDev1"}, nbins);
y2 = y2 / sum(y2);

% Display the Histograms
figure(4); hold on; grid on; 
plot(x1,y1,'r',x2,y2,'b');
xlabel("std values");
ylabel("Relative frequency");
title("Histogram of standard Deviation for EOG E1-M2")
legend("REM", "Non-REM");

% ----------------------------------------------------------------

% Histogram of standard deviation (EOG E2-M2) for REM sleep
[y1, x1] = hist(features{mask1 & mask3, "stdDev2"}, nbins);
y1 = y1 / sum(y1);

% Histogram of standard deviation (EOG E2-M2) for NREM sleep
[y2, x2] = hist(features{mask2 & mask3, "stdDev2"}, nbins);
y2 = y2 / sum(y2);

% Display the Histograms
figure(5); hold on; grid on; 
plot(x1,y1,'r',x2,y2,'b');
xlabel("std values");
ylabel("Relative frequency");
title("Histogram of standard Deviation for EOG E2-M2")
legend("REM", "Non-REM");
