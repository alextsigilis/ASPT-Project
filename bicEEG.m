% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 17th, 2022
% -------------------------------------------------------------------
%
% Function Description:
% This function estimates the bicoherence index of 
% EEG recordings in 30sec windows. By locating the 
% peaks of the bicoherence, we can detect quadratic
% phase coupling between different frequencies. The
% bicoherence is essentially a normalized bispectrum.
% Many methods have been proposed to normalize the 
% bispectrum. This function uses the one proposed by
% Hagihira 2001 and Hayashi 2007. 
% -------------------------------------------------------------------
%
% Arguments List: (X, K, fs, fc, channel)
%
% X: (table) the EEG recordings and sleep stage
% Annotations. You should use loadEDF to obtain this 
% table.
%
% K: (integer) number of partitions for estimating 
% the bicoherence in a 30sec epoch. The bispectrum is 
% estimated for every partition. The bicoherence is 
% obtained by normalizing the bispectra across the entire
% 30sec epoch.
%
% fs: (float) sampling frequency of EEG recordings
% (Usually 256Hz).
%
% fc: (float) The upper limit of the frequency axis.
% Since it is rarely necessary to estimate the
% bispectrum in the for every pair of frequencies
% -fs/2 < f1,f2 < +fs/2, we can greatly reduce the 
% execution time and memory requirements by truncating
% the frequency axis. Always make sure that fs and 
% fc satisfy the Nyquist criterion. That is: 2*fc < fs
%
% channel: (integer) An integer between 1 and 4 used to 
% select one out of the 4 available EEG channels.
% -------------------------------------------------------------------
%
% Return Variables: (bic, freq)
%
% bic: (table) A table with two columns. The first column
% stores the bicoherence estimations for every 30sec epoch.
% The second column stores the sleep stage Annotations 
%
% freq: (array of floats) A 1D array which can be used as a 
% 
% ===================================================================

function [bic, freq] = bicEEG(X, fs, fc, K, channel)

    % ---------------------------------------------------------------
    % Parameter Checks
    % ---------------------------------------------------------------  
    
    if K  <= 0 error("K must be positive "); end
    if fs <= 0 error("fs must be positive"); end
    if fc <= 0 error("fc must be positive"); end
    if channel < 1 || channel > 4 error("invalid EEG channel"); end 
    
    N = size(X,1);
    if N <= 0 error("X is empty"); end
    
    L = numel(cell2mat(X{1,channel}));
    if L <= 0 error("EEG records are empty"); end
    
    M = floor(L/K);
    if M < 5 error("low frequency resolution"); end
    if 2*M*fc >= (M-2)*fs error("Nyquist criterion is violated"); end
    if fc*M <= fs error("fc is almost zero"); end
    
    % ---------------------------------------------------------------
    % Create a table to store the bicoherence and sleep stage labels
    % ---------------------------------------------------------------
    
    names = [X.Properties.VariableNames{channel}, "Annotations"];
    types = ["cell", "string"];
    
    bic = table(                    ...
        'Size',          [N 2],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);
    
    bic.Annotations = X.Annotations;
    
    % ---------------------------------------------------------------
    % Frequency axis
    % ---------------------------------------------------------------
    
    if (rem(M,2) == 0)
        freq = [-M/2:(M/2-1)];
    else
        freq = [-(M-1)/2:(M-1)/2];
    end

    % ---------------------------------------------------------------
    % Estimation of Bispectra
    % ---------------------------------------------------------------

    % idx: array of indices for obtaining a truncated FFT
    % len: length of truncated FFT
    idx = 1:M;
    idx = idx(-fc * M <= freq * fs & freq * fs <= fc * M);
    len = numel(idx);

    % tri: array of indices for accumulating triple products
    tri = hankel([1:len],[len,1:len-1] );

    % Bispectrum estimations for every 30sec record
    for i = 1:1:N
        % x:   (1D array) a 30sec EEG recording
        % ind: (1D array) array of indices for partitioning x
        % win: (1D array) time-domain Hanning window
        ind = [1:M]';
        x   = cell2mat(X{i,channel});
        win = hanning(M);

        % B,F3,F12: (2D matrices) temporary variables
        % for accumulating triple products
        B   = zeros(len,len);
        F12 = zeros(len,len);
        F3  = zeros(len,len);

        % Calculate partial bispectra for every partition
        for j = 1:1:K
            % y: (1D array) a partition of x
            % Y: (1D array) the FFT of y
            y = x(ind); y = y - mean(y);                                         
            Y = fft(y .* win);
            
            % Discard any frequencies above +fc and below -fc 
            % by truncating the result of the previous FFT
            Y = fftshift(Y);
            Y = Y(idx);                                                
            Y = ifftshift(Y);
            CY = conj(Y);

            % Estimate triple products
            temp1 = reshape(CY(tri),len,len);
            temp2 = Y * Y.';
            
            % Accumulate triple products
            B   = B   + temp1 .* temp2;
            F3  = F3  + abs(temp1) .^ 2;
            F12 = F12 + abs(temp2) .^ 2;
            
            % update the indices of ind
            ind = ind + M;                                    
        end

        % Calculate the bicoherence indices
        b = (abs(B) .^ 2) ./ (F12 .* F3);

        % Shift the elements of the bicoherence matrix
        b = fftshift(b);
        
        % Save the results and proceed to the next record
        bic{i,1} = {b};
    end

    % ---------------------------------------------------------------
    % Truncate the frequency axis and rescale according to fs 
    % ---------------------------------------------------------------
    freq = freq(idx);
    freq = freq * fs / M;
end