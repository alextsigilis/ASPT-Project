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

% Initialize empty tables to store bicoherence features
sz = [0 6];
types = ["double" "double" "double" "double" "double" "string"];
names = ["ent1" "ent2" "ent3" "avg" "peak" "Annotations"];
delta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
theta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
alpha = table('Size',sz,'VariableTypes',types,'VariableNames',names);
beta  = table('Size',sz,'VariableTypes',types,'VariableNames',names);

% Extract bicoherence features from every patient
for i = start:stop
    txt = sprintf("SN%03d.edf",i);
    if ~isfile(txt) continue; end

    fprintf("Loading EDF files for patient %d ... ",i);
    Z = loadEDF(i);
    fprintf("OK\n");

    fprintf("Estimating bicoherence matrix ... ");
    [X,f] = bicEEG(Z,K,fs,fc,channel);
    fprintf("OK\n");

    fprintf("Extracting bicoherence features ... ");
    [D, U, A, B] = bicoherFeatures(X,f);
    delta = [delta; D];
    theta = [theta; U];
    alpha = [alpha; A];
    beta  = [beta;  B];
    fprintf("OK\n");

    fprintf("\n");
end

idx = 1;

% --------------- Histograms for delta-wave features -----------------

W  = delta.Annotations == "Sleep stage W";
R  = delta.Annotations == "Sleep stage R";
N1 = delta.Annotations == "Sleep stage N1";
N2 = delta.Annotations == "Sleep stage N2";
N3 = delta.Annotations == "Sleep stage N3";

K = size(delta,2) - 1;

for k = 1:K
    z1 = delta{W,  k};
    z2 = delta{R,  k};
    z3 = delta{N1, k};
    z4 = delta{N2, k};
    z5 = delta{N3, k};

    [y1, x1] = hist(z1,nbins); y1 = y1 / numel(z1);
    [y2, x2] = hist(z2,nbins); y2 = y2 / numel(z2);
    [y3, x3] = hist(z3,nbins); y3 = y3 / numel(z3);
    [y4, x4] = hist(z4,nbins); y4 = y4 / numel(z4);
    [y5, x5] = hist(z5,nbins); y5 = y5 / numel(z5);

    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(delta.Properties.VariableNames{k}); 
    ylabel('Probability');
    title('Delta waves');
    legend("W", "R", "N1", "N2", "N3");
end

% --------------- Histograms for theta-wave features -----------------

W  = theta.Annotations == "Sleep stage W";
R  = theta.Annotations == "Sleep stage R";
N1 = theta.Annotations == "Sleep stage N1";
N2 = theta.Annotations == "Sleep stage N2";
N3 = theta.Annotations == "Sleep stage N3";

K = size(theta,2) - 1;

for k = 1:K
    z1 = theta{W,  k};
    z2 = theta{R,  k};
    z3 = theta{N1, k};
    z4 = theta{N2, k};
    z5 = theta{N3, k};

    [y1, x1] = hist(z1,nbins); y1 = y1 / numel(z1);
    [y2, x2] = hist(z2,nbins); y2 = y2 / numel(z2);
    [y3, x3] = hist(z3,nbins); y3 = y3 / numel(z3);
    [y4, x4] = hist(z4,nbins); y4 = y4 / numel(z4);
    [y5, x5] = hist(z5,nbins); y5 = y5 / numel(z5);

    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(theta.Properties.VariableNames{k}); 
    ylabel('Probability');
    title('theta waves');
    legend("W", "R", "N1", "N2", "N3");
end

% --------------- Histograms for alpha-wave features -----------------

W  = alpha.Annotations == "Sleep stage W";
R  = alpha.Annotations == "Sleep stage R";
N1 = alpha.Annotations == "Sleep stage N1";
N2 = alpha.Annotations == "Sleep stage N2";
N3 = alpha.Annotations == "Sleep stage N3";

K = size(alpha,2) - 1;

for k = 1:K
    z1 = alpha{W,  k};
    z2 = alpha{R,  k};
    z3 = alpha{N1, k};
    z4 = alpha{N2, k};
    z5 = alpha{N3, k};

    [y1, x1] = hist(z1,nbins); y1 = y1 / numel(z1);
    [y2, x2] = hist(z2,nbins); y2 = y2 / numel(z2);
    [y3, x3] = hist(z3,nbins); y3 = y3 / numel(z3);
    [y4, x4] = hist(z4,nbins); y4 = y4 / numel(z4);
    [y5, x5] = hist(z5,nbins); y5 = y5 / numel(z5);

    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(alpha.Properties.VariableNames{k}); 
    ylabel('Probability');
    title('alpha waves');
    legend("W", "R", "N1", "N2", "N3");
end

% --------------- Histograms for beta-wave features -----------------

W  = beta.Annotations == "Sleep stage W";
R  = beta.Annotations == "Sleep stage R";
N1 = beta.Annotations == "Sleep stage N1";
N2 = beta.Annotations == "Sleep stage N2";
N3 = beta.Annotations == "Sleep stage N3";

K = size(beta,2) - 1;

for k = 1:K
    z1 = beta{W,  k};
    z2 = beta{R,  k};
    z3 = beta{N1, k};
    z4 = beta{N2, k};
    z5 = beta{N3, k};

    [y1, x1] = hist(z1,nbins); y1 = y1 / numel(z1);
    [y2, x2] = hist(z2,nbins); y2 = y2 / numel(z2);
    [y3, x3] = hist(z3,nbins); y3 = y3 / numel(z3);
    [y4, x4] = hist(z4,nbins); y4 = y4 / numel(z4);
    [y5, x5] = hist(z5,nbins); y5 = y5 / numel(z5);

    figure(idx); idx = idx + 1;
    plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5); grid on;
    xlabel(beta.Properties.VariableNames{k}); 
    ylabel('Probability');
    title('beta waves');
    legend("W", "R", "N1", "N2", "N3");
end