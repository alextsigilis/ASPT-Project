% ======================================================
% Demo Script (3)
% ------------------------------------------------------
%
% Author: Christodoulos Michaelides, July 2022
% ------------------------------------------------------
%
% Objectives:
% ------------------------------------------------------
% 1) Use Multi-Resolution-Analysis on the entire 
% recording of an EEG/EOG/ECG channel
%
% 2) Choose the most appropriate frequency scales
% to reconstruct an approximation of the original
% recording
% 
% 3) Estimate 2nd, 3rd and 4th order moments
% for every frequency scale independently. 
% Use a small time window of 1-3 seconds, when 
% estimating variance, skewness and kurtosis in 
% order to obtain localised estimates of the above
% quantities.
% 
% 4) See if there is any correlation between 
% variance/kurtosis/skewness and transitional phenomena
% (K-complex, arousal periods) on the reconstructed EEG
% by visual inspection of the waveforms.
% ======================================================

% reset your workspace and clear terminal
clear all; close all; clc;

% ======================================================
% 1) Parameters of MRA Decomposition and Reconstruction
% (Choose any value that you want to experiment with)
% ======================================================

% name of input file
input_file = "SN001.edf";            

% name of annotations file
annot_file = "SN001_sleepscoring.edf";

% index of selected channel
% index: signal label
% 1:     EEG F4-M1
% 2:     EEG C4-M1
% 3:     EEG O2-M1
% 4:     EEG C3-M2
% 5:     EMG chin
% 6:     EOG E1-M2
% 7:     EOG E2-M2
% 8:     ECG
channel = 1;

% window lenth for estimating higher order statistics
% (you should probably leave this as it is)
w = duration("00:00:30");

% wavelet for MRA
wavelet = "db2";

% number of levels on MRA binary tree
levels = 6;                          

% logical array for reconstruction 
% -> must be levels+1 elements long,
%    otherwise the script won't run properly
% -> frequency boundaries are calculated
%    assuming a sampling rate of 256Hz
levelForReconstruction = [
    false,  ...     % scale 1 (64-128 Hz)
    false,  ...     % scale 2 (32-64 Hz)  
    true,   ...     % scale 3 (16-32 Hz)
    true,   ...     % scale 4 (8-16 Hz)
    true,   ...     % scale 5 (4-8 Hz)
    true,   ...     % scale 6 (2-4 Hz)
    true    ...     % scale 7 (0-2 Hz)
];                   
                                     
% array of valid sleep stages
% (you should probably leave this 
% as it is, unless you want to add
% new sleep stages)
stages = cellstr(       ...
    ["Sleep stage W",   ...
     "Sleep stage N1",  ...
     "Sleep stage N2",  ...
     "Sleep stage N3",  ...
     "Sleep stage R"]   ...
);

% ------------------------------------------------------
% Do not change anything below that point, 
% unless you know exactly what you are doing.
% ------------------------------------------------------

% number of frequency scales
num_of_scales = levels + 1;            

% validity checks on selected frequency scales
if length(levelForReconstruction) ~= num_of_scales
    fprintf('Error:\n');
    fprintf('Invalid Selection of Frequency Scales\n');
    fprintf('Abort ...\n');
    return;
end

if all(levelForReconstruction == false)
    fprintf('Error:\n');
    fprintf('Select at least one frequency scale\n');
    fprintf('Abort ...\n');
    return;
end

% ======================================================
% 2) Read data from the EDF file
% ======================================================

% Progress Status
fprintf("Loading input file ...  "); tic;

% Read the entire EDF file
X = edfread(input_file);

% choose the appropriate channel
sig = X{:,channel};
sig = cell2mat(sig);     

% Read the annotations from the 
% sleep scoring EDF file
[~,labels] = edfread(annot_file);

% Read file metadata
info = edfinfo(input_file);

% N:  number of data records per recording
% d:  duration of every data record in seconds
% n:  samples per data record
% fs: sampling frequency in Hertz
N = info.NumDataRecords;
d = seconds(info.DataRecordDuration);
n = info.NumSamples(channel);
fs = n / d;

% construct a time axis for
% the EEG/EOG/ECG channel
time = linspace(0, N*d, N*d*fs);  

% Remove Unnecessary Columns and Rows from Annotations
labels = removevars(labels, "Duration");
rows = ~ismember(categorical(labels.Annotations), categorical(stages)); 
labels(rows,:) = [];

% Join the timetable of the EEG/EOG/ECG recordings
% with the timetable of the sleep labels
Y = synchronize(X, labels, X.("Record Time"), "previous");
labels = Y.("Annotations");
labels = categorical(labels);
labels = reordercats(labels, stages);
labels = renamecats(labels, stages, ["W" "N1" "N2" "N3" "R"]);

% delete unused variables
clear X Y;

% Progress Status
fprintf("Done\n\n"); toc; fprintf("\n");

% ======================================================
% 3) Perform MRA decomposition and reconstruction
% ======================================================

% Progress Status
fprintf("Performing MRA ... "); tic;

% Perform the decomposition using modwt
wt = modwt(sig,wavelet,levels);

% Construct MRA matrix using modwtmra
mra = modwtmra(wt,wavelet);

% Sum down the rows of the selected multiresolution signals
sig1 = sum(mra(levelForReconstruction,:),1);

% Delete unused variables
clear wt;  

% Progress Status
fprintf("Done\n\n"); toc; fprintf("\n"); 

% ======================================================
% 4) Plot of Original vs Reconstructed Waveform
% ======================================================            

% Reconstructed vs Original Signal
figure(1); 

plt = plot(time,sig, time,sig1);
plt(1).LineWidth = 0.5; 
plt(2).LineWidth = 2.0;

legend('Original Signal', 'Reconstructed Signal');

xlabel("time in seconds");
ylabel("Amplitude in microVolts");
title("Original vs Reconstruction");

% ======================================================
% Estimate variance/skewness/kurtosis
% in a sliding window.
% ======================================================

% length of sliding window in samples
w = seconds(w);
l = w * fs;

% number of segments created by sliding window
K = N * d * fs / l; K = floor(K);

% frequency boundaries for every analysis scale
freq = (fs/2) * 2.^-[0:levels];
freq = [freq 0];

% Plots of selected frequency scales
% and sliding window statistics
for i = 1:1:num_of_scales
    % ignore discarded frequency scales
    if levelForReconstruction(i) == false
        continue;
    end

    % Progress Status
    fprintf("Estimating statistics for scale %d ... ", i);
    
    % f1: lower bound of frequency scale in Hertz
    % f2: upper bound of frequency scale in Hertz
    f1 = freq(i+1);
    f2 = freq(i+0);

    % estimate higher order moments
    x = mra(i,1:K*l);
    x = reshape(x, K, l);
    m   = sum(x,2) / l;            % sliding average 
    var = sum((x-m).^2, 2) / l;    % sliding variance
    std = sqrt(var);               % sliding standard deviation
    skw = sum((x-m).^3, 2) / l;    % sliding 3rd moment
    krt = sum((x-m).^4, 2) / l;    % sliding 4th moment
    skw = skw ./ (std.^3);         % sliding skewness
    krt = krt ./ (std.^4);         % sliding kurtosis

    % Scatter plot of high order statistics
    figure(2*i);
    plot3(std, skw, krt); grid on;
    xlabel('Standard Deviation');
    ylabel('Skewness');
    zlabel('Kurtosis');
    title(sprintf('Scatter Plot: %.2fHz - %.2fHz', f1, f2));

    % Plots of higher order statistics
    % with respect to time
    figure(2*i+1); t = tiledlayout(5,1);

    % upsample to match the length
    % of the frequency scale
    std = interp1(1:1:K, std, (1:1:K*l)/l);
    skw = interp1(1:1:K, skw, (1:1:K*l)/l);
    krt = interp1(1:1:K, krt, (1:1:K*l)/l);

    % subplot of histogram
    x0 = nexttile;
    plot(labels);
    xlabel('time in seconds');
    ylabel('Sleep stage');
    title('Hypnogram');

    % subplot of current frequency scale
    x1 = nexttile;
    plot(x1, time(1:K*l), mra(i,1:K*l));
    xlabel('time in seconds'); 
    ylabel(sprintf('Scale %d', i));
    title(sprintf('%.2fHz-%.2fHz', f1, f2));
    
    % subplot of sliding variance
    x2 = nexttile;
    plot(x2, time(1:K*l), std);
    xlabel('time in seconds');
    ylabel('Std. Deviation');

    % subplot of sliding skewness
    x3 = nexttile;             
    plot(x3, time(1:K*l), skw);
    xlabel('time in seconds');
    ylabel('Skewness');

    % subplot of sliding kurtosis
    x4 = nexttile;                     
    plot(x4, time(1:K*l), krt); 
    xlabel('time in seconds');
    ylabel('Kurtosis');

    % Link time axes
    linkaxes([x1 x2 x3 x4], 'x');

    % Progress Status
    fprintf("Done\n\n");
end