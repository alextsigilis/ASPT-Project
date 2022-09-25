% ================================================================
% Author: Christodoulos Michaelides
% Date: September 24th, 2022
% ----------------------------------------------------------------
%
% Script Description
%
% 1) Choose a single patient and EEG/ECG channel from the entire
% dataset
%
% 2) Estimate the cepstral coefficients for every 30sec epoch
%
% 3) Display plots of the estimated cepstra
%
% 4) Save the plots as .fig or .png files for later use
% ================================================================

% ================================================================
% Reset your Matlab workspace
% ================================================================
clear all; close all; clc;

% ================================================================
% Script parameters. Choose any value you want to experiment with
% ================================================================

% patient_ID: (int) selected patient
% channel: (int or string) selected EEG/ECG recording
patient_ID = 3;
channel = "EEGC3_M2";

% dt: (int or float) duration of time window for cepstrum
% fs: (int or float) sampling frequency of EEG/ECG recording
dt = 5.0;
fs = 256;

% save folder for cepstrum plots
% Adjust this variable to be compatible with your filesystem
path = sprintf("C:\\Users\\User\\Desktop\\cepFigures");

% ------------- Do not change anything below that point ------------- %

% ================================================================
% Create empty folders to save the cepstrum plots
% ================================================================

% W:  (string) sub-folder for Sleep stage W
% R:  (string) sub-folder for Sleep stage R
% N1: (string) sub-folder for Sleep stage N1
% N2: (string) sub-folder for Sleep stage N2
% N3: (string) sub-folder for Sleep stage N3
W  = sprintf("%s\\Sleep stage W",  path);
R  = sprintf("%s\\Sleep stage R",  path);
N1 = sprintf("%s\\Sleep stage N1", path);
N2 = sprintf("%s\\Sleep stage N2", path);
N3 = sprintf("%s\\Sleep stage N3", path);

% Delete previous save-folders
if isfolder(path)
    rmdir(path,'s');
end

% Create new save-folders
mkdir(path);
mkdir(W);
mkdir(R);
mkdir(N1);
mkdir(N2);
mkdir(N3);

% ================================================================
% Load EEG/ECG recordings and extract cepstral coefficients
% ================================================================

% Load ECG/EEG recordings
fprintf("Loading EEG/ECG recordings ... ");
Z = loadEDF(patient_ID);
fprintf("Done\n");

% Estimate cepstral coefficients
fprintf("Estimating cepstral coefficients ... ");
[C, t] = cepEEG(Z, fs, dt, channel);
fprintf("Done\n");

% ================================================================
% Generate cepstrum plots and save them 
% ================================================================

K = size(C,1);

% Plot every vector of cepstral coefficients
fprintf("Generating plots of cepstral coefficients ...\n");
pause(3); clc;

for k = 1:1:K
    % Progress bar
    fprintf("Progress: %d/%d\n",k,K);

    % Plot of cepstral coefficients
    f = figure(1);
    plot(t, cell2mat(C{k,"ceps"})); grid on;

    % Crop the horizontal and vertical axes
    xlim([0 dt/2]);
    ylim([-5 +5]);

    % Add axis-labels and title
    xlabel("quefrency in sec");
    ylabel("cepstral coefficients");
    title(sprintf("Plot %d, %s",k, C{k,"Annotations"}));

    % Save the cepstrum plot
    savepath = sprintf("%s\\%s\\%d.png", path, C{k,"Annotations"}, k);
    saveas(f, savepath);

    % Pause the execution for 0.1 seconds.
    % This is required because of a bug in
    % Matlab R2022b which occasionally leads
    % to random crashes.
    % Increase the pause-duration, if those
    % random crashes persist.
    pause(0.1); clc;
end

