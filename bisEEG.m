% ===================================================================
% Author: Christodoulos Michaelides
% Date: August 17th, 2022
% -------------------------------------------------------------------
%
% Function Description:
% Use this function to estimate the bispectra of EEG recordings
% in 30sec windows (by applying the direct estimation method). This 
% implementation is heavily influenced by the Higher Order Spectral
% Analysis Toolbox (HOSA Toolbox). The main benefit of our 
% implementation is that it allows the user to estimate the 
% bispectrum in a pre-defined frequency region. By carefully 
% choosing the boundaries of this region, the execution time and 
% memory requirements can be dramatically reduced without losing any
% significant amount of information from the original bispectrum. 
% -------------------------------------------------------------------
%
% Arguments List: (X, K, fs, fc, channel)
%
% X: (table) the EEG recordings and sleep stage Annotations. Use
% loadEDF to generate this table. You can process multiple EEG tables
% at once by using loadEDF on many patients and (vertically) 
% concatenating their EEG tables (assuming that your system can 
% handle the extra memory requirements).
% 
% K: (integer) The bispectrum is obtained by splitting every 30sec 
% EEG record into K non-overlapping segments. A partial bispectrum is 
% estimated for every segment. The final estimation is obtained by 
% averaging the partial bispectra. Bear in mind that using a large 
% value for K results in poor frequency resolution. Make sure that 
% every segment has at least 5 EEG samples, otherwise an error 
% message is displayed and the estimation is terminated. A good value
% for K is K = 4
%
% fs: (float) The sampling frequency of the EEG recordings in Hertz 
% (typically 256Hz)
% 
% fc: (float) The upper limit of the bispectrum frequency axes 
% in Hertz. Since it is not necessary to estimate the bispectrum 
% for every pair of frequencies:  -fs/2 < f1,f2 < fs/2, 
% we can significantly reduce the execution time and memory 
% requirements by truncating the bispectrum. 
% A good choice for fc is fc = 6Hz. Always make sure that fc and
% fs do not violate the Nyquist criterion. That is 2*fc < fs
%
% channel: (integer) An integer between 1 and 4 used to select 
% one out of the 4 EEG channels for every patient.
% -------------------------------------------------------------------
%
% Return Variables: (bis, freq)
%
% bis: (table) A table with two columns. The first column contains
% the bispectra for every 30sec EEG recording. The second column
% contains the corresponding sleep stage Annotation.
%
% freq: (array of floats) A 1D array which can be used as a frequency
% axis when plotting the bispectra. The upper and lower bound is 
% approximately -fc and fc respectively.
% ===================================================================

function [bis, freq] = bisEEG (X, K, fs, fc, channel)
    % ---------------------------------------------------------------
    % Parameter Checks
    % ---------------------------------------------------------------

    if K  <= 0 error("K must be positive "); end                           
    if fs <= 0 error("fs must be positive"); end                           
    if fc <= 0 error("fc must be positive"); end                                       
    if channel < 1 || channel > 4 error("invalid EEG channel"); end        

    % N: (integer) number of 30sec EEG recordings
    N = size(X,1);
    if N <= 0 error("X is empty"); end 

    % L: (integer) number of samples per 30sec recording
    L = numel(cell2mat(X{1,channel}));
    if L <= 0 error("EEG records are empty"); end

    % M: (integer) number of samples per subset after partitioning
    M = floor(L/K);                                                                                                                                                            
    if M < 5 error("low frequency resolution"); end
    if 2*M*fc >= (M-2)*fs error("Nyquist criterion is violated"); end  
    if fc * M <= fs error("fc is almost zero"); end

    % ---------------------------------------------------------------
    % Construction of frequency domain window function
    % ---------------------------------------------------------------
    
    % W: (integer) Size of Rao-Gabr window (must be an odd number)
    W = 5;
    
    % temporary variables for constructing the Rao-Gabr window
    Q = fix(M/W);
    u = (-(W-1)/2):(+(W-1)/2);
    
    % opwind: Rao-Gabr window for frequency domain smoothing
    opwind = ones(W,1) * (u .^2);
    opwind = opwind + opwind' + u' * u;
    opwind = 1 - (2*Q/M)^2 * opwind;
    
    % hex: (boolean matrix) hexagonal boolean mask
    hex = ones(W,1) * u;
    hex = abs(hex) + abs(hex') + abs(hex+hex');
    hex = (hex < W);
    
    % Apply the hexagonal mask
    % Rescale the remaining values appropriately
    opwind = opwind .* hex;
    opwind = opwind ./ sum(opwind(:));

    % ---------------------------------------------------------------
    % Create a table to store bispectra and sleep stage labels
    % ---------------------------------------------------------------
    
    names = [X.Properties.VariableNames{channel}, "Annotations"];
    types = ["cell", "string"];
    
    bis = table(                    ...
        'Size',          [N 2],     ...
        'VariableTypes', types,     ...
        'VariableNames', names);
    
    bis.Annotations = X.Annotations;

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

    % Hexagonal mask to remove artifacts outside
    % the hexagonal symmetry regions.
    lo = -(len-1)/2; hi = +(len-1)/2; u = lo:1:hi;
    hex = ones(len,1) * u;
    hex = abs(hex) + abs(hex') + abs(hex+hex');
    hex = (hex < len);

    % tri: array of indices for accumulating triple products
    tri = hankel([1:len],[len,1:len-1] );

    % Bispectrum estimations for every 30sec record
    for i = 1:1:N
        % x: a 30sec EEG recording
        % B: 2D matrix for storing the bispectrum of x
        % ind: array of indices for extracting partitions from x
        B  = zeros(len,len);                                              
        ind = [1:M]';
        x = cell2mat(X{i,channel});

        % Calculate partial bispectra for every partition
        for j = 1:1:K
            % y: (1D array) a partition of x
            % Y: (1D array) the FFT of y
            y = x(ind);                                             
            Y = fft(y-mean(y)) / M;
            
            % Discard any frequencies above +fc and below -fc 
            % by truncating the result of the previous FFT
            Y = fftshift(Y);
            Y = Y(idx);                                                
            Y = ifftshift(Y);

            % Accumulate triple products
            CY = conj(Y);
            B  = B + (Y * Y.') .* reshape(CY(tri), len, len);    
            
            % update the indices of ind
            ind = ind + M;                                    
        end

        % Shift the elements of the bispectrum matrix
        B = fftshift(B) / K;

        % Frequency domain smoothing with the Rao-Gabr window
        B = conv2(B,opwind,'same');    
        
        % Save the results and proceed to the next record
        bis{i,1} = {B .* hex};
    end

    % ---------------------------------------------------------------
    % Truncate the frequency axis and rescale it according to fs 
    % ---------------------------------------------------------------
    freq = freq(idx);
    freq = freq * fs / M;
return