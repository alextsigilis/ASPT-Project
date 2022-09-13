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
useDWT = false;     % Extract features from DWT 
                    % coefficients or time-domain

% ==========================================================
% Extract features from EOG/EMG recordings
% ==========================================================

sz = [0 6];
types = ["double" "double" "double" "double" "double" "string"];
names = ["xcorr" "ECB1" "ECB2" "stdEMG" "VppEMG" "Annotations"];
features = table('Size',sz,'VariableTypes',types,'VariableNames',names);

for i = id1:1:id2
    edf = sprintf("SN%03d.edf",i);
    mat = sprintf("%03d.mat",i);

    if (~isfile(edf)) && (~isfile(mat)) continue; end

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

isAwake = features{:,"Annotations"} == "Sleep stage W";

% ----------------------------------------------------------------

% Histogram of cross-correlation for REM sleep
data = features{isREM, "xcorr"};
[y1, x1] = hist(data, nbins);
y1 = y1 / sum(y1);

% Histogram of cross-correlation for NREM sleep 
data = features{isNREM, "xcorr"};
[y2, x2] = hist(data, nbins);
y2 = y2 / sum(y2);

% Histogram of cross-correlation for Sleep stage W
data = features{isAwake, "xcorr"};
[y3, x3] = hist(data, nbins);
y3 = y3 / sum(y3);

% Display the histograms
figure(1); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b',x3,y3,'k');
xlabel("EOG cross-correlation coefficients");
ylabel("probability");
title("Histogram of EOG cross-correlation");
legend("REM","Non-REM","Awake");

% ----------------------------------------------------------------

% Histogram of EMG variance for Sleep stage R
data = features{isREM,"stdEMG"}; data = log10(1 + data);
[y1, x1] = hist(data, nbins);
y1 = y1 / sum(y1);

% Histogram of EMG variance for Sleep stages N1, N2 and N3
data = features{isNREM,"stdEMG"}; data = log10(1 + data);
[y2, x2] = hist(data, nbins);
y2 = y2 / sum(y2);

% Histogram of EMG variance for Sleep stage W
data = features{isAwake,"stdEMG"}; data = log10(1 + data);
[y3, x3] = hist(data, nbins);
y3 = y3 / sum(y3);

% Display the histograms
figure(2); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b',x3,y3,'k');
xlabel("Standard deviation of chin-EMG (log scale)");
ylabel("probability");
title("Histogram of EMG standard deviation");
legend("REM","Non-REM","Awake");

% ----------------------------------------------------------------

% Histogram of ECB1 for Sleep stage R
data = features{isREM,"ECB1"}; data = log10(data);
[y1, x1] = hist(data, nbins);
y1 = y1 / sum(y1);

% Histogram of ECB1 for Sleep stages N1, N2 and N3
data = features{isNREM,"ECB1"}; data = log10(data);
[y2, x2] = hist(data, nbins);
y2 = y2 / sum(y2);

% Histogram of ECB1 for Sleep stage W
data = features{isAwake,"ECB1"}; data = log10(data);
[y3, x3] = hist(data, nbins);
y3 = y3 / sum(y3);

% Display the histograms
figure(3); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b',x3,y3,'k');
xlabel("ECB1, log scale");
ylabel("probability");
title("Histogram of ECB for 1st EOG channel");
legend("REM", "NREM", "Awake");

% ----------------------------------------------------------------

% Histogram of ECB2 for Sleep stage R
data = features{isREM,"ECB2"}; data = log10(data);
[y1, x1] = hist(data, nbins);
y1 = y1 / sum(y1);

% Histogram of ECB2 for Sleep stages N1, N2 and N3
data = features{isNREM, "ECB2"}; data = log10(data);
[y2, x2] = hist(data, nbins);
y2 = y2 / sum(y2);

% Histogram of ECB2 for Sleep stage W
data = features{isAwake,"ECB2"}; data = log10(data);
[y3, x3] = hist(data, nbins);
y3 = y3 / sum(y3);

% Display the histograms
figure(4); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b',x3,y3,'k');
xlabel("ECB2, log scale");
ylabel("probability");
title("Histogram of ECB for 2nd EOG channel");
legend("REM", "NREM", "Awake");

% ----------------------------------------------------------------

% Histogram of EMG peak-to-peak voltage for Sleep stage R
data = features{isREM,"VppEMG"}; data = log10(data);
[y1, x1] = hist(data, nbins);
y1 = y1 / sum(y1);

% Histogram of EMG peak-to-peak voltage for Sleep stages N1, N2 and N3
data = features{isNREM,"VppEMG"}; data = log10(data);
[y2, x2] = hist(data, nbins);
y2 = y2 / sum(y2);

% Histogram of EMG peak-to-peak voltage for Sleep stage W
data = features{isAwake,"VppEMG"}; data = log10(data);
[y3, x3] = hist(data, nbins);
y3 = y3 / sum(y3);

% Display the histograms
figure(5); hold on; grid on;
plot(x1,y1,'r',x2,y2,'b',x3,y3,'k');
xlabel("Peak-to-peak voltage of chin-EMG (log scale)");
ylabel("probability");
title("Histogram of EMG peak-to-peak voltage");
legend("REM", "NREM", "Awake");

