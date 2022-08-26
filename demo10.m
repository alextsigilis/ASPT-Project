% Clear all workspace variables,
% Close any open windows and
% clear the command line.
clear all; close all; clc;

% ------------------------ Script Parameters ------------------------

% Hyperparameters for 
% estimating the bicoherence
K = 32;                         % Number of segments
fs = 256;                       % Sampling frequency
fc = 32;                        % upper bound on frequency axis

% Select patients and an EEG
% channels from the dataset
channel = 1;                    % EEG channel
start = 1;                      % first patient
stop = 50;                      % last patient

% Plot settings
nbins = 50;                     % Number of histogram bins

% ---------------- Do not change anything below ---------------------

% Initialize an empty table to store bicoherence features
sz = [0 4];
types = ["double" "double" "double" "string"];
names = ["ent1" "ent2" "ent3" "Annotations"];
Y = table('Size',sz,'VariableTypes',types,'VariableNames',names);

% Extract bicoherence features from every patient
for i = start:stop
    txt = sprintf("SN%03d.edf",i);
    if ~isfile(txt) continue; end

    fprintf("Loading EDF files for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("OK\n");

    fprintf("Estimating bicoherence matrix ... ");
    [X,~] = bicEEG(Z,K,fs,fc,channel);
    fprintf("OK\n");

    fprintf("Extracting bicoherence features ... ");
    Y = [Y; bicoherFeatures(X)];
    fprintf("OK\n");

    fprintf("\n");
end

% ---------- Plot histograms of the bicoherence features ------------

% Boolean masks to distinguish between the five sleep stages
W  = Y.Annotations == "Sleep stage W";
R  = Y.Annotations == "Sleep stage R";
N1 = Y.Annotations == "Sleep stage N1";
N2 = Y.Annotations == "Sleep stage N2";
N3 = Y.Annotations == "Sleep stage N3";

% histograms of bicoherence entropy
[y1,x1] = hist(Y{W,"ent1"},nbins);   y1 = y1 / sum(y1);
[y2,x2] = hist(Y{R,"ent1"},nbins);   y2 = y2 / sum(y2);
[y3,x3] = hist(Y{N1,"ent1"},nbins);  y3 = y3 / sum(y3);
[y4,x4] = hist(Y{N2,"ent1"},nbins);  y4 = y4 / sum(y4);
[y5,x5] = hist(Y{N3,"ent1"},nbins);  y5 = y5 / sum(y5);

figure(1); grid on;
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence entropy");
ylabel("probability");
legend("W", "R", "N1", "N2", "N3");

% histogram of bicoherence squared-entropy
[y1,x1] = hist(Y{W,"ent2"},nbins);   y1 = y1 / sum(y1);
[y2,x2] = hist(Y{R,"ent2"},nbins);   y2 = y2 / sum(y2);
[y3,x3] = hist(Y{N1,"ent2"},nbins);  y3 = y3 / sum(y3);
[y4,x4] = hist(Y{N2,"ent2"},nbins);  y4 = y4 / sum(y4);
[y5,x5] = hist(Y{N3,"ent2"},nbins);  y5 = y5 / sum(y5);

figure(2); grid on;
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence squared entropy");
ylabel("probability");
legend("W", "R", "N1", "N2", "N3");

% histograms of bicoherence cubic-entropy
[y1,x1] = hist(Y{W,"ent3"},nbins);   y1 = y1 / sum(y1);
[y2,x2] = hist(Y{R,"ent3"},nbins);   y2 = y2 / sum(y2);
[y3,x3] = hist(Y{N1,"ent3"},nbins);  y3 = y3 / sum(y3);
[y4,x4] = hist(Y{N2,"ent3"},nbins);  y4 = y4 / sum(y4);
[y5,x5] = hist(Y{N3,"ent3"},nbins);  y5 = y5 / sum(y5);

figure(3); grid on;
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence cubed entropy");
ylabel("probability");
legend("W", "R", "N1", "N2", "N3");