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
% method: (string) choose between real and complex cepstrum 
%   => "rceps" for real cepstrum coefficients
%   => "cceps" for complex cepstrum coefficients
dt = 5;
fs = 256;
method = "rceps";

% path: (string) save folder for cepstrum plots. Adjust 
% this variable to be compatible with your filesystem
path = sprintf("C:\\Users\\User\\Desktop\\cepFigures");

% format: (string) file extension for saving the cepstrum plots.
% You can choose between "fig" and "png". The first option ("fig") is 
% a proprietary file-format provided by Matlab which allows us to 
% open the plots later by using the builtin openfig() MATLAB command.
% The second option ("png") is an lossless-compression image format 
% suitable for graphics objects.
format = "fig";

% limits: (2x1 array of floats) Maximum and minimum values for the 
% y-axis of the cepstrum plots.  You can experiment until you find
% a suitable pair of values.
limits = [0 +1e-2]; 

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
[C, t] = cepEEG(Z, fs, dt, channel, method);
fprintf("Done\n");

% ================================================================
% Generate cepstrum plots and save them 
% ================================================================

% Notch filter for removing power-line interference
f0 = 50; Q  = 64;
w0 = 2*f0 / fs; bw = w0 / Q;
[num, den] = iirnotch(w0,bw);

% Plot every vector of cepstral coefficients
fprintf("Generating plots of cepstral coefficients ...\n");
pause(3); clc;

K = size(C,1);

for k = 1:1:K
    % Progress bar
    fprintf("Progress: %d/%d\n",k,K);

    % New display window
    f = figure(1);

    % Plot of cepstral coefficients
    subplot(2,1,1);
    plot(t, cell2mat(C{k,"ceps"}).^2); grid on;

    % Crop the horizontal and vertical axes
    xlim([0 dt/2]);
    ylim(limits);

    % Add axis-labels and title
    xlabel("quefrency in seconds");
    ylabel("squared cepstral coefficients");
    title(sprintf("Plot %d, %s",k, C{k,"Annotations"}));

    % Plot of original EEG/ECG recording 
    subplot(2,1,2);
    sig = cell2mat(Z{k,channel});
    sig = filtfilt(num,den,double(sig));
    tau = linspace(0,30,numel(sig));
    plot(tau, sig); grid on;

    % Add axis-labels and title
    xlabel("time in seconds");
    ylabel("signal in uVolts");

    % Save the cepstrum plot
    savepath = sprintf("%s\\%s\\%d.%s", path, C{k,2}, k, format);
    
    if format == "fig" || format == "png"
        saveas(f, savepath);
    else
        error("Invalid file format");
    end

    % Pause the execution for 0.1 seconds.
    % This is required because of a bug in
    % Matlab R2022b which occasionally leads
    % to random crashes.
    % Increase the pause-duration, if those
    % random crashes persist.
    pause(0.1); clc;
end