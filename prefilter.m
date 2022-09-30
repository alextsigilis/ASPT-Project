% ================================================================
% Authors: Stergios Gregoriou, Christodoulos Michaelides
% Date: September 28th, 2022
% ----------------------------------------------------------------
%
% Function Description:
% The process of recording EEG signals is susceptible to 
% various noise sources such as powerline interference and
% electrode drift.
% 
% Powerline interference manifests itself as a large spike in the 
% power-spectrum of the recorded signal located at 50Hz or 60Hz 
% depending on the frequency of the AC.
%
% Electrode drift refers to low frequency components in the
% EEG recordings (typically below 0.1Hz). The cause of these
% artifacts is usually eye-blinking and bad electrode contacts.
% 
% The purpose of this function is to suppress those artifacts
% as much as possible.
% ----------------------------------------------------------------
% 
% Function Arguments:
%
% X: (table) The polysomnographic recordings of a patient. You 
% should use loadEDF() to obtain this table.
%
% fs: (int) The sampling frequency of the recordings expressed
% in Hertz.
% 
% channel: (string) The name of the channel which we want to 
% filter.
% ---------------------------------------------------------------- 
%
% Return Variables:
%
% Y: (table) A table with two columns. The first column contains
% 30sec epochs of the filtered signal. The second column contains
% the Sleep stage Annotations.
% ================================================================

function Y = prefilter(X, fs, channel)
    % H: an elliptical filter designed to remove
    % powerline interference at 50Hz and electrode
    % drift below 0.2Hz
    H = designfilt(                         ...
        'bandpassiir',                      ...
        'StopbandFrequency1',   0.1,        ...
        'PassbandFrequency1',   0.4,        ...
        'PassbandFrequency2',   33.5,       ...
        'StopbandFrequency2',   50,         ...
        'StopbandAttenuation1', 20,         ...
        'StopbandAttenuation2', 20,         ...
        'PassbandRipple',       0.01,       ...
        'SampleRate',           fs,         ...
        'DesignMethod',         'ellip');

    % N: number of 30sec epochs
    % M: number of samples per 30sec epoch
    N = size(X,1); sz = [N 2];
    M = numel(cell2mat(X{1,channel}));

    % Y: a table which holds the 
    % output of the elliptic filter.
    names = [channel "Annotations"];
    types = ["cell" "string"];

    Y = table(                      ...
        'Size',          sz,        ...
        'VariableTypes', types,     ...
        'VariableNames', names);

    % Copy sleep stage Annotations
    Y.Annotations = X.Annotations;
    
    % Extract the entire EEG signal
    y = cell2mat(X{:,channel}); y = y(:);

    % Zero padding
    y = [zeros(fs, 1); y; zeros(fs, 1)];

    % Apply elliptical filter
    y = filtfilt(H,double(y));
    y = single(y);

    % Remove zero padding;
    y = y((fs+1):1:end-fs);

    % Split to 30sec epochs
    y = reshape(y, [M N]);

    % Save the output of the filter
    for n = 1:1:N
        Y{n,channel} = {y(:,n)};
    end
end