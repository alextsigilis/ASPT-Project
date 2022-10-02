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
sz = [0 6]; K = 5;
types = ["single" "single" "single" "single" "single" "string"];
names = ["mean" "std" "skw" "krt" "zcr" "Annotations"];

Delta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
Theta = table('Size',sz,'VariableTypes',types,'VariableNames',names);
Alpha = table('Size',sz,'VariableTypes',types,'VariableNames',names);
Beta  = table('Size',sz,'VariableTypes',types,'VariableNames',names);

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
    Z = prefilter(Z,256,channel);

    % Extract cepstral features
    fprintf("Extracting cepstral features ... \n");
    [delta, theta, alpha, beta] = rcepFeatures(Z,channel);

    Delta = [Delta; delta];
    Theta = [Theta; theta];
    Alpha = [Alpha; alpha];
    Beta =  [Beta;  beta ];

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

% binary masks for delta wave cepstral features
W  = Delta.Annotations == "Sleep stage W";
R  = Delta.Annotations == "Sleep stage R";
N1 = Delta.Annotations == "Sleep stage N1";
N2 = Delta.Annotations == "Sleep stage N2";
N3 = Delta.Annotations == "Sleep stage N3";

% Histograms of delta-wave cepstral features
for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Delta{W,i};
    x2 = Delta{R,i};
    x3 = Delta{N1,i};
    x4 = Delta{N2,i};
    x5 = Delta{N3,i};

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
    xlabel(names(i));
    ylabel("Probability density function");
    title("Delta waves / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for theta wave cepstral features
W  = Theta.Annotations == "Sleep stage W";
R  = Theta.Annotations == "Sleep stage R";
N1 = Theta.Annotations == "Sleep stage N1";
N2 = Theta.Annotations == "Sleep stage N2";
N3 = Theta.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Theta{W,i};
    x2 = Theta{R,i};
    x3 = Theta{N1,i};
    x4 = Theta{N2,i};
    x5 = Theta{N3,i};

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
    xlabel(names(i));
    ylabel("Probability density function");
    title("Theta waves / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for alpha wave cepstral features
W  = Alpha.Annotations == "Sleep stage W";
R  = Alpha.Annotations == "Sleep stage R";
N1 = Alpha.Annotations == "Sleep stage N1";
N2 = Alpha.Annotations == "Sleep stage N2";
N3 = Alpha.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Alpha{W,i};
    x2 = Alpha{R,i};
    x3 = Alpha{N1,i};
    x4 = Alpha{N2,i};
    x5 = Alpha{N3,i};

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
    xlabel(names(i));
    ylabel("Probability density function");
    title("Alpha waves / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end

% Binary masks for beta wave cepstral features
W  = Beta.Annotations == "Sleep stage W";
R  = Beta.Annotations == "Sleep stage R";
N1 = Beta.Annotations == "Sleep stage N1";
N2 = Beta.Annotations == "Sleep stage N2";
N3 = Beta.Annotations == "Sleep stage N3";

for i = 1:1:K
    % Split feature table based on sleep stage annotations
    x1 = Beta{W,i};
    x2 = Beta{R,i};
    x3 = Beta{N1,i};
    x4 = Beta{N2,i};
    x5 = Beta{N3,i};

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
    xlabel(names(i));
    ylabel("Probability density function");
    title("Beta waves / Cepstral features");
    saveas(f,sprintf("%s\\%d.fig", path, idx)); idx = idx + 1;
end