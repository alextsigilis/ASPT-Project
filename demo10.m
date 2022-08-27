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
nbins = 100;                    % Number of histogram bins
norm = "pdf";                   % Normalization type for histograms:
                                % cdf => cumulative distribution 
                                % pdf => probability density function

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
z1 = Y{W, "ent1"}; z1 = log10(1+z1);
z2 = Y{R, "ent1"}; z2 = log10(1+z2);
z3 = Y{N1,"ent1"}; z3 = log10(1+z3);
z4 = Y{N2,"ent1"}; z4 = log10(1+z4);
z5 = Y{N3,"ent1"}; z5 = log10(1+z5);

figure(1); hold on; grid on;
[y1,x1] = hist(z1,nbins); y1 = y1 / numel(z1);
[y2,x2] = hist(z2,nbins); y2 = y2 / numel(z2);
[y3,x3] = hist(z3,nbins); y3 = y3 / numel(z3);
[y4,x4] = hist(z4,nbins); y4 = y4 / numel(z4);
[y5,x5] = hist(z5,nbins); y5 = y5 / numel(z5);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence entropy (log scale)");
ylabel("probability");
legend("W","R","N1","N2","N3");

fprintf("Statistics for bicoherence entropy:\n");
fprintf("-------------------------------------------\n");
fprintf("stage W  => mean: %.4f std: %.4f\n", mean(z1), std(z1));
fprintf("stage R  => mean: %.4f std: %.4f\n", mean(z2), std(z2));
fprintf("stage N1 => mean: %.4f std: %.4f\n", mean(z3), std(z3));
fprintf("stage N2 => mean: %.4f std: %.4f\n", mean(z4), std(z4));
fprintf("stage N3 => mean: %.4f std: %.4f\n", mean(z5), std(z5));
fprintf("\n");

% histogram of bicoherence squared-entropy
z1 = Y{W, "ent2"}; z1 = log10(1+z1);
z2 = Y{R, "ent2"}; z2 = log10(1+z2);
z3 = Y{N1,"ent2"}; z3 = log10(1+z3);
z4 = Y{N2,"ent2"}; z4 = log10(1+z4);
z5 = Y{N3,"ent2"}; z5 = log10(1+z5);

figure(2); hold on; grid on;
[y1,x1] = hist(z1,nbins); y1 = y1 / numel(z1);
[y2,x2] = hist(z2,nbins); y2 = y2 / numel(z2);
[y3,x3] = hist(z3,nbins); y3 = y3 / numel(z3);
[y4,x4] = hist(z4,nbins); y4 = y4 / numel(z4);
[y5,x5] = hist(z5,nbins); y5 = y5 / numel(z5);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence squared-entropy (log scale)");
ylabel("probability");
legend("W","R","N1","N2","N3");

fprintf("Statistics for bicoherence squared-entropy:\n");
fprintf("-------------------------------------------\n");
fprintf("stage W  => mean: %.4f std: %.4f\n", mean(z1), std(z1));
fprintf("stage R  => mean: %.4f std: %.4f\n", mean(z2), std(z2));
fprintf("stage N1 => mean: %.4f std: %.4f\n", mean(z3), std(z3));
fprintf("stage N2 => mean: %.4f std: %.4f\n", mean(z4), std(z4));
fprintf("stage N3 => mean: %.4f std: %.4f\n", mean(z5), std(z5));
fprintf("\n");

% histograms of bicoherence cubic-entropy
z1 = Y{W, "ent3"}; z1 = log10(1+z1);
z2 = Y{R, "ent3"}; z2 = log10(1+z2);
z3 = Y{N1,"ent3"}; z3 = log10(1+z3);
z4 = Y{N2,"ent3"}; z4 = log10(1+z4);
z5 = Y{N3,"ent3"}; z5 = log10(1+z5);

figure(3); hold on; grid on;
[y1,x1] = hist(z1,nbins); y1 = y1 / numel(z1);
[y2,x2] = hist(z2,nbins); y2 = y2 / numel(z2);
[y3,x3] = hist(z3,nbins); y3 = y3 / numel(z3);
[y4,x4] = hist(z4,nbins); y4 = y4 / numel(z4);
[y5,x5] = hist(z5,nbins); y5 = y5 / numel(z5);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5);
xlabel("bicoherence cubic-entropy (log scale)");
ylabel("probability");
legend("W","R","N1","N2","N3");

fprintf("Statistics for bicoherence cubic-entropy:\n");
fprintf("-------------------------------------------\n");
fprintf("stage W  => mean: %.4f std: %.4f\n", mean(z1), std(z1));
fprintf("stage R  => mean: %.4f std: %.4f\n", mean(z2), std(z2));
fprintf("stage N1 => mean: %.4f std: %.4f\n", mean(z3), std(z3));
fprintf("stage N2 => mean: %.4f std: %.4f\n", mean(z4), std(z4));
fprintf("stage N3 => mean: %.4f std: %.4f\n", mean(z5), std(z5));
fprintf("\n");