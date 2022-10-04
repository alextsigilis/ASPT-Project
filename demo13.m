% ===================================================================
% Authors:  Chrysa Doulou,
%           Demetrios Orphanos,
%           Christodoulos Michaelides
%   
% Date:     October 1st, 2022     
% -------------------------------------------------------------------
%   
% Script Description:
% ---------------------
% 1) Choose a number of patients from the dataset and an EEG channel.
% 2) Split the EEG recordings into delta, theta, alpha and beta wave
%    components.
% 3) Estimate the cepstal sequences of every wave type.
% 4) Extract features from the estimated cepstra
% 5) Plot histograms to see if there is any relationship between 
%    those features and the sleep stage annotations.
% ===================================================================

% Reset your Matlab Workspace
clear all; close all; clc;

% ===================================================================
% Script parameters. (Choose any values you want to experiment with)
% ===================================================================

% Choose a number of patients and an EEG channel
first = 1;
last  = 25;
channel = "EEGC3_M2";
fs = 256;

% Number of bins for histograms
nbins = 40;

% Directory for saving generated histograms
path = sprintf("C:\\Users\\USER\\Desktop\\cepHist");

% ------------ Do not change anything below that point --------------

% ===================================================================
% Extraction of Cepstral Features
% ===================================================================

% Initialize empty tables to store cepstral features and sleep 
% stage annotations.
sz1 = [0 4];
sz2 = [0 4];
types1 = ["single" "single" "single" "string"];
types2 = ["single" "single" "single" "string"];
names1 = ["std" "skw" "krt" "Annotations"];
names2 = ["mean" "std" "zcr" "Annotations"];

K = 3;

Alpha1 = table('Size',sz1,'VariableTypes',types1,'VariableNames',names1);
Alpha2 = table('Size',sz2,'VariableTypes',types2,'VariableNames',names2);
Beta1  = table('Size',sz1,'VariableTypes',types1,'VariableNames',names1);
Beta2  = table('Size',sz2,'VariableTypes',types2,'VariableNames',names2);

for n = first:1:last
    
    edfFile = sprintf("SN%03d.edf",n);
    matFile = sprintf("%03d.mat",n);

    if (~isfile(edfFile)) && (~isfile(matFile))
        continue;
    end

    % Load EEG recordings
    fprintf("Patient %d\n",n);
    Z = loadEDF(n);

    % Prefiltering
    fprintf("Prefiltering ...\n");
    Z = prefilter(Z,fs);

    % Extract cepstral features
    fprintf("Extracting cepstral features ... \n");
    [alpha1, beta1, alpha2, beta2] = rcepFeatures(Z,channel);

    Alpha1 = [Alpha1; alpha1];
    Alpha2 = [Alpha2; alpha2];
    Beta1  = [Beta1;  beta1];
    Beta2  = [Beta2;  beta2];

    fprintf("\n\n");
end

% ===================================================================
% Histograms of extracted features
% ===================================================================

% Delete previous save folder
if isfolder(path)
    rmdir(path,'s');
end

% Make a new save folder
mkdir(path);

% Index for histogram plots
idx = 1;

% binary masks for cepstral features
W  = Alpha1.Annotations == "Sleep stage W";
R  = Alpha1.Annotations == "Sleep stage R";
N1 = Alpha1.Annotations == "Sleep stage N1";
N2 = Alpha1.Annotations == "Sleep stage N2";
N3 = Alpha1.Annotations == "Sleep stage N3";

% Histograms of cepstral features
for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Alpha1{W,i};
    x2 = Alpha1{R,i};
    x3 = Alpha1{N1,i};
    x4 = Alpha1{N2,i};
    x5 = Alpha1{N3,i};

    % Construct histograms for every class-sleep stage annotation
    [h1, t1] = hist(x1,nbins);
    [h2, t2] = hist(x2,nbins);
    [h3, t3] = hist(x3,nbins);
    [h4, t4] = hist(x4,nbins); 
    [h5, t5] = hist(x5,nbins);

    % Normalize histograms to obtain pdf-estimations
    h1 = (h1 * nbins) / (numel(x1) * range(x1));
    h2 = (h2 * nbins) / (numel(x2) * range(x2));
    h3 = (h3 * nbins) / (numel(x3) * range(x3));
    h4 = (h4 * nbins) / (numel(x4) * range(x4));
    h5 = (h5 * nbins) / (numel(x5) * range(x5));

    f = figure;
    plot(t1,h1,t2,h2,t3,h3,t4,h4,t5,h5); grid on;
    legend("W", "R", "N1", "N2", "N3");
    xlabel(names1(i));
    ylabel("Probability density function");
    title("Alpha waves (1) / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for alpha wave cepstral features
W  = Alpha2.Annotations == "Sleep stage W";
R  = Alpha2.Annotations == "Sleep stage R";
N1 = Alpha2.Annotations == "Sleep stage N1";
N2 = Alpha2.Annotations == "Sleep stage N2";
N3 = Alpha2.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Alpha2{W,i};
    x2 = Alpha2{R,i};
    x3 = Alpha2{N1,i};
    x4 = Alpha2{N2,i};
    x5 = Alpha2{N3,i};

    % Construct histograms for every class-sleep stage annotation
    [h1, t1] = hist(x1,nbins);
    [h2, t2] = hist(x2,nbins);
    [h3, t3] = hist(x3,nbins);
    [h4, t4] = hist(x4,nbins); 
    [h5, t5] = hist(x5,nbins);

    % Normalize histograms to obtain pdf-estimations
    h1 = (h1 * nbins) / (numel(x1) * range(x1));
    h2 = (h2 * nbins) / (numel(x2) * range(x2));
    h3 = (h3 * nbins) / (numel(x3) * range(x3));
    h4 = (h4 * nbins) / (numel(x4) * range(x4));
    h5 = (h5 * nbins) / (numel(x5) * range(x5));

    f = figure;
    plot(t1,h1,t2,h2,t3,h3,t4,h4,t5,h5); grid on;
    legend("W", "R", "N1", "N2", "N3");
    xlabel(names2(i));
    ylabel("Probability density function");
    title("Alpha waves (2) / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for beta wave cepstral features
W  = Beta1.Annotations == "Sleep stage W";
R  = Beta1.Annotations == "Sleep stage R";
N1 = Beta1.Annotations == "Sleep stage N1";
N2 = Beta1.Annotations == "Sleep stage N2";
N3 = Beta1.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Beta1{W,i};
    x2 = Beta1{R,i};
    x3 = Beta1{N1,i};
    x4 = Beta1{N2,i};
    x5 = Beta1{N3,i};

    % Construct histograms for every class-sleep stage annotation
    [h1, t1] = hist(x1,nbins);
    [h2, t2] = hist(x2,nbins);
    [h3, t3] = hist(x3,nbins);
    [h4, t4] = hist(x4,nbins); 
    [h5, t5] = hist(x5,nbins);

    % Normalize histograms to obtain pdf-estimations
    h1 = (h1 * nbins) / (numel(x1) * range(x1));
    h2 = (h2 * nbins) / (numel(x2) * range(x2));
    h3 = (h3 * nbins) / (numel(x3) * range(x3));
    h4 = (h4 * nbins) / (numel(x4) * range(x4));
    h5 = (h5 * nbins) / (numel(x5) * range(x5));

    f = figure;
    plot(t1,h1,t2,h2,t3,h3,t4,h4,t5,h5); grid on;
    legend("W", "R", "N1", "N2", "N3");
    xlabel(names1(i));
    ylabel("Probability density function");
    title("Beta waves (1) / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for beta wave cepstral features
W  = Beta2.Annotations == "Sleep stage W";
R  = Beta2.Annotations == "Sleep stage R";
N1 = Beta2.Annotations == "Sleep stage N1";
N2 = Beta2.Annotations == "Sleep stage N2";
N3 = Beta2.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Beta2{W,i};
    x2 = Beta2{R,i};
    x3 = Beta2{N1,i};
    x4 = Beta2{N2,i};
    x5 = Beta2{N3,i};

    % Construct histograms for every class/Sleep stage Annotation
    [h1, t1] = hist(x1,nbins);
    [h2, t2] = hist(x2,nbins);
    [h3, t3] = hist(x3,nbins);
    [h4, t4] = hist(x4,nbins);
    [h5, t5] = hist(x5,nbins);

    % Normalize histograms to obtain pdf-estimations
    h1 = (h1 * nbins) / (numel(x1) * range(x1));
    h2 = (h2 * nbins) / (numel(x2) * range(x2));
    h3 = (h3 * nbins) / (numel(x3) * range(x3));
    h4 = (h4 * nbins) / (numel(x4) * range(x4));
    h5 = (h5 * nbins) / (numel(x5) * range(x5));

    f = figure;
    plot(t1,h1,t2,h2,t3,h3,t4,h4,t5,h5); grid on;
    legend("W", "R", "N1", "N2", "N3");
    xlabel(names2(i));
    ylabel("Probability density function");
    title("Beta waves (2) / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end