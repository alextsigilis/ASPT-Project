% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 5th, 2022
% -------------------------------------------------------------------
%
% Demo Script 6:
%
% 1) Choose one out of the 150 patients in the dataset
%
% 2) Select one out of the 4 available EEG channels
%
% 3) Choose five 30sec recordings out of the selected EEG channel
%    Ideally there should be one recording for every Sleep Stage.
%    That is:
%    one 30sec recording for Sleep stage W
%    one 30sec recording for Sleep stage N1
%    one 30sec recording for Sleep stage N2
%    one 30sec recording for Sleep stage N3
%    one 30sec recording for Sleep stage R
%
% 4) Estimate the Bispectra of those EEG recordings.
%
% 5) Plot the magnitudes of the estimated bispectra in a contour
%    plot
%
% 6) Examine the results by visual inspection of the plots.
% ------------------------------------------------------------------

clear all; close all; clc;

% ==================================================================
% 1) Script Parameters: Choose any value you want to experiment with
% ==================================================================

% Patient parameters
patient_ID  = 1;    % patient index   (1 ... 151)
EEG_channel = 2;    % EEG channel index (1 ... 4)

% Sleep stage parameters
idx1 = 1;           % time index for Sleep stage W
idx2 = 13;          % time index for Sleep stage N1
idx3 = 67;          % time index for Sleep stage N2
idx4 = 647;         % time index for Sleep stage N3
idx5 = 532;         % time index for Sleep stage R

% Parameters for bispectrum estimation
K  = 3;
fs = 256;       % EEG sampling frequency
fc = 4;         % bispectrum truncation frequency

% Parameters for bispectrum contour plots
levels = 16;    % Number of contour lines for every plot

% ==================================================================
% 2) Estimate the bispectrum
% ==================================================================

fprintf("Loading EEG recordings ... ");
X = loadEDF(patient_ID);
fprintf("Done\n");

fprintf("Estimating Bispectra ... ");
[bis, freq] = bisEEG(X, K, fs, fc, EEG_channel);
fprintf("Done\n");

% ==================================================================
% 3) Plot the selected bispectra
% ==================================================================

B1 = cell2mat(bis{idx1,1});    % Extract the bispectrum for Sleep stage W
B2 = cell2mat(bis{idx2,1});    % Extract the bispectrum for Sleep stage N1
B3 = cell2mat(bis{idx3,1});    % Extract the bispectrum for Sleep stage N2
B4 = cell2mat(bis{idx4,1});    % Extract the bispectrum for Sleep stage N3
B5 = cell2mat(bis{idx5,1});    % Extract the bispectrum for Sleep stage R

figure(1); grid on;
contourf(freq, freq, abs(B1), levels, 'LineColor', 'none');
xlabel('f_1'); ylabel('f_2'); title('Sleep stage W');

figure(2); grid on; 
contourf(freq, freq, abs(B2), levels, 'LineColor', 'none');
xlabel('f_1'); ylabel('f_2'); title('Sleep stage N1');

figure(3); grid on;
contourf(freq, freq, abs(B3), levels, 'LineColor', 'none');
xlabel('f_1'); ylabel('f_2'); title('Sleep stage N2');

figure(4); grid on;
contourf(freq, freq, abs(B4), levels, 'LineColor', 'none');
xlabel('f_1'); ylabel('f_2'); title('Sleep stage N3');

figure(5); grid on;
contourf(freq, freq, abs(B5), levels,'LineColor','none');
xlabel('f_1'); ylabel('f_2'); title('Sleep stage R');