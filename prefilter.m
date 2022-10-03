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
% ---------------------------------------------------------------- 
%
% Return Variables:
%
% Y: (table) A table with N+1 columns. The first N column contain
% 30sec epochs of the filtered signals. The last column contains
% the Sleep stage Annotations.
% ================================================================

function Y = prefilter(X, fs)
    % H: an elliptic filter designed to remove
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

    % N:  number of PSG channels
    % K: number of 30sec epochs per channel
    % dt: duration of epochs in seconds
    % M:  number of samples per 30sec epoch
    N  = size(X,2) - 1;
    K  = size(X,1);
    dt = 30;
    M  = fs * dt;

    for n = 1:1:N
        % Extract an entire PSG channel
        x = cell2mat(X{:,n}); 
        x = x(:);

        % Zero padding
        x = [zeros(fs,1); x; zeros(fs,1)];
        
        % Apply elliptic filter
        x = filtfilt(H,double(x)); x = single(x);

        % Remove zero padding
        x = x((fs+1):1:(end-fs));

        % Split channel to 30sec epochs again
        x = reshape(x, [M K]);

        % Store the output of the filter
        for k = 1:1:K
            X{k,n} = {x(:,n)};
        end
    end

    Y = X;
end