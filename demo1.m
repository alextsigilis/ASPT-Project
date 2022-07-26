% =========================================================
% Demo Script (1)
%
% Objectives:
% ------------------------------------
% 1) Reading Data from EDF/EDF+ files
% 2) Using Matlab's timetable data type
% =========================================================

% =========================================================
% 1) File I/O for the EDF/EDF+ file format. 
% =========================================================

% Read recordings from an EDF file
X = edfread("SN001.edf");

% Use edfinfo to get additional
% information about the EDF file
info = edfinfo("SN001.edf");

% Get the name of every recording in the EDF file.
labels = info.SignalLabels;

% Get the number of recordings in the EDF file.
N = info.NumSignals;

% Get the number of DataRecords per signal
L = info.NumDataRecords;

% Get the number of samples per DataRecord
M = info.NumSamples;

% =========================================================
% 2) Using MATLAB's timetable data structure.
% =========================================================

% Column Indexing with curly braces and label-names
% (This is the recommended way of indexing columns)
% ---------------------------------------------------------
F4M1 = X{:,"EEGF4_M1"};           % 1st EEG channel
C4M1 = X{:,"EEGC4_M1"};           % 2nd EEG channel
O2M1 = X{:,"EEGO2_M1"};           % 3rd EEG channel
C3M2 = X{:,"EEGC3_M2"};           % 4th EEG channel
EMG  = X{:,"EMGChin"};            % chin EMG channel
E1M2 = X{:,"EOGE1_M2"};           % 1st EOG channel
E2M2 = X{:,"EOGE2_M2"};           % 2nd EOG channel
ECG  = X{:,"ECG"};                % ECG channel

% Column Indexing with curly braces and integers
% (works just fine but it's not the recommended way)
% ---------------------------------------------------------
% F4M1 = X{:,1};
% C4M1 = X{:,2};
% O2M1 = X{:,3};
% C3M2 = X{:,4};
% EMG  = X{:,5};
% E1M2 = X{:,6};
% E2M2 = X{:,7};
% ECG  = X{:,8};

% Column Indexing with Parentheses and label-names
% ---------------------------------------------------------
% F4M1 = X.("EEGF4_M1"); 
% C4M1 = X.("EEGC4_M1");
% O2M1 = X.("EEGO2_M1");
% C3M2 = X.("EEGC3_M2");
% EMG  = X.("EMGChin");
% E1M2 = X.("EOGE1_M2");
% E2M2 = X.("EOGE2_M2");
% ECG  = X.("ECG");

% Row Indexing with curly braces and integers
% ---------------------------------------------------------
sig1 = X{1:30,"EEGF4_M1"};  % first 30 seconds of the F4-M1

% Row Indexing with curly braces and duration()/timerange()
% ---------------------------------------------------------
t0 = duration(0,0,0);           % start at t0 = 0sec
dt = duration(0,0,30);          % move 30 seconds forward
interval  = range(t0,t0+dt);    % [0sec, 30sec)  

sig2 = X{interval,"EEGF4_M1"};  % first 30 seconds of F4-M1

% Converting cell arrays into ordinary 1D-arrays
% ---------------------------------------------------------
sig3 = cell2mat(sig2);