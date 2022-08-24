% ==========================================================
% Author: Christodoulos Michaelides
% Date: August 22nd, 2022
% ----------------------------------------------------------
%
% Demo Script 8:
%
% 1) Choose a number of patients from the dataset
% 2) Extract features from the EOG and EMG recordings of 
%    those patients
% 3) Plot histograms of those features for various sleep 
%    stages.
% ==========================================================

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
useDWT = true;      % Extract features from DWT 
                    % coefficients or time-domain

% ==========================================================
% Extract features from EOG recordings
% ==========================================================

sz = [0 5];
types = ["double" "double" "double" "double" "string"];
names = ["xcorr" "ECB1" "ECB2" "stdEMG" "Annotations"];
features = table('Size',sz,'VariableTypes',types,'VariableNames',names);

for i = id1:1:id2
    if ~isfile(sprintf("SN%03d.edf",i)) continue; end

    fprintf("Loading EOG/EMG recordings for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("Done\n");

    fprintf("Extracting features from EOG/EMG recordings ... ");
    features = [features; features_EOG_EMG(Z,useDWT)];
    fprintf("Done\n");
end

% ==========================================================
% Plot histograms/scatterplots of features
% ==========================================================

% isREM:    (1D array) boolean array for REM epochs
% isNREM:   (1D array) boolean array for NREM epochs
% isAsleep: (1D array) boolean array for REM or NREM epochs 

isREM = features{:,"Annotations"} == "Sleep stage R";

isNREM =                                                ... 
    (features{:,"Annotations"} == "Sleep stage N1") |   ...
    (features{:,"Annotations"} == "Sleep stage N2") |   ...
    (features{:,"Annotations"} == "Sleep stage N3");

isAsleep = isREM | isNREM;

% ----------------------------------------------------------------

% Histogram of cross-correlation for REM sleep
[y1, x1] = hist(features{isREM, "xcorr"}, nbins);
y1 = y1 / sum(y1);

% Histogram of cross-correlation for NREM sleep 
[y2, x2] = hist(features{isNREM, "xcorr"}, nbins);
y2 = y2 / sum(y2);

% Display the histograms
figure(1); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b');
xlabel("EOG cross-correlation coefficients");
ylabel("probability");
title("Histogram of EOG cross-correlation");
legend("REM","Non-REM");

% ----------------------------------------------------------------

% Histogram of EMG variance for Sleep stage W
[y1, x1] = hist(features{~isAsleep, "stdEMG"}, 10*nbins);
y1 = y1 / sum(y1);

% Histogram of EMG variance for Sleep stages R, N1, N2 and N3
[y2, x2] = hist(features{ isAsleep, "stdEMG"}, 10*nbins);
y2 = y2 / sum(y2);

% Display the histograms
figure(2); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b');
xlabel("Standard deviation of chin-EMG");
ylabel("probability");
title("Histogram of EMG standard deviation");
legend("Sleep stage W", "Sleep stages R,N1,N2 and N3");

% ----------------------------------------------------------------

% Histogram of ECB1 for Sleep stage W
[y1, x1] = hist(features{isREM, "ECB1"}, 10*nbins);
y1 = y1 / sum(y1);

% Histogram of ECB1 for Sleep stages R, N1, N2 and N3
[y2, x2] = hist(features{isNREM, "ECB1"}, 10*nbins);
y2 = y2 / sum(y2);

% Display the histograms
figure(3); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b');
xlabel("ECB1");
ylabel("probability");
title("Histogram of ECB for 1st EOG channel");
legend("REM", "NREM");

% ----------------------------------------------------------------

% Histogram of ECB2 for Sleep stage W
[y1, x1] = hist(features{isREM, "ECB2"}, 10*nbins);
y1 = y1 / sum(y1);

% Histogram of ECB2 for Sleep stages R, N1, N2 and N3
[y2, x2] = hist(features{isNREM, "ECB2"}, 10*nbins);
y2 = y2 / sum(y2);

% Display the histograms
figure(4); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b');
xlabel("ECB2");
ylabel("probability");
title("Histogram of ECB for 2nd EOG channel");
legend("REM", "NREM");