% Orfanos Dimitrios, 9579

% Function Description:
% This function can be used to read EEG recordings and
% sleep scorings from EDF files and convert them into
% MATLAB timetables.
% ---------------------------------------------------------
%
% Function Arguments:
% idx: (integer) the index of the patient's heart rate whose 
% recordings we want to read. There are 151 patients in the dataset.
% Every patient is given an index ranging from 1 to 154.
% Certain indices are missing. For example if we are
% interested in reading the recordings of patient 4, then
% we would use the following command:
% >>> Z = loadEDF_HeartRateECG(4);
%
% The Heart Rate exported with the EDFbrowser. 
%
% ---------------------------------------------------------
%
% Return Variables:
% Z: (timetable) A timetable which contains the heart rate
% and the sleep stage labels for the selected patient.
%
% The "Onset" column is the time-axis of the timetable.
% Every element on this column corresponds to the onset
% of a 30second labeled segment in the recordings.
% The "Annotations" column contains the sleep stage labels
% for every segment.
% =========================================================

function Z = loadEDF_HeartRateECG(idx)
% =====================================================
% 1) Read data from the EDF files
% =====================================================

% input_file: (string) the name of the file
% which contains the EEG recordings
% annot_file: (string) the name of the file
% which contains the EEG annotations
input_file = sprintf("SN%03d_HR.edf",idx);
annot_file = sprintf("SN%03d_sleepscoring.edf",idx);

% Make sure that the input files exist
if ~isfile(input_file) || ~isfile(annot_file)
    error("Error: Input files not found\n\n");
end

% X:    (timetable) the EEG/EOG/EMG/ECG recordings
% y:    (timetable) the sleep stage annotations
% info: metadata of the input file
[X,~] = edfread(input_file);
[~,y] = edfread(annot_file);
info  = edfinfo(input_file);

% N: (integer) number of data records per recording
% M: (integer) number of channels per patient
% d: (integer) duration of every record in seconds (1 sec)
% n: (integer) number of samples per record (256 samples)
N = info.NumDataRecords;
M = info.NumSignals;
d = seconds(info.DataRecordDuration);
n = info.NumSamples;
channel_names = X.Properties.VariableNames;

% =====================================================
% 2) Removal of unnecessary of sleep-stage labels
% =====================================================

% A cell array of every valid sleep
% stage that we expect to find in the
% annotations file.
% ---------------------------------------
% "Sleep stage W":     awake
% "Sleep stage N1":    stage 1 NREM sleep
% "Sleep stage N2":    stage 2 NREM sleep
% "Sleep stage N3":    stage 3 NREM sleep
% "Sleep stage R"      REM sleep
valid_stages = cellstr([ ...
    "Sleep stage W",     ...
    "Sleep stage N1",    ...
    "Sleep stage N2",    ...
    "Sleep stage N3",    ...
    "Sleep stage R"]);

% w: the duration of every labeled segment in seconds
w = duration("00:00:30");

% Find rows which contain invalid labels.
% Certain labels contain information about
% the light conditions inside the recording room.
% Those labels should be discarded
rows = ~ismember(                ...
    categorical(y.Annotations),  ...
    categorical(valid_stages)    ...
    );

% remove invalid rows
y(rows,:) = [];

% Make sure that the remaining labels
% correspond to 30 seconds intervals
if (any(y{:,"Duration"} ~= w))
    error("Error: Invalid Annotations Detected\n\n");
end

% remove unused columns. Since every labeled
% segment corresponds to a 30 second interval,
% we can safely remove the "Duration" column.
y = removevars(y,"Duration");

% =====================================================
% 3) Synchronize the EEG and Annotations timetables
% The Annotations timetable (y) uses a 30sec step for
% its time axis. The EEG timetable (X) uses a 1sec
% for its time axis. Therefore, it is necessary to
% synchronize the time axes of the timetables before
% merging them.
% =====================================================

% K: (integer) number of 30sec long segments per channel
% Z: (array) synchronized timetable
[K,~] = size(y);
Z = cell(K,M);

for i = 1:1:K
    % t0: onset of current label in seconds
    % [t0,t1]: time interval of current label
    t0 = y.Onset(i);
    t1 = t0 + w;
    dt = timerange(t0,t1);

    for j = 1:1:M
        x = cell2mat(X{dt,j});
        Z(i,j) = num2cell(x,1);
    end
end

% convert cell-array to timetable and adjust variable names
Z = cell2table(Z);
Z = renamevars(Z,Z.Properties.VariableNames,channel_names);
Z = addvars(Z,y.Annotations,'NewVariableNames','Annotations');
Z = addvars(Z,y.Onset,'NewVariableNames','Onset');
Z = table2timetable(Z);
end