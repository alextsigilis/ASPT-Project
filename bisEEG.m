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

    if K  <= 0 error("K must be positive");  end                           
    if fs <= 0 error("fs must be positive"); end                           
    if fc <= 0 error("fc must be positive"); end                                       
    if channel < 1 || channel > 4 error("invalid EEG channel"); end        

    nrecs = size(X,1);
    if nrecs <= 0 error("X is empty"); end 
    ly = numel(cell2mat(X{1,channel}));
    if ly <= 0 error("EEG records are empty"); end

    nsamp = floor(ly/K);                                                   
    nfft = nsamp;                                                          
    nadvance = nsamp;                                                      

    if nsamp < 5 error("low frequency resolution"); end
    if 2*nfft*fc >= (nfft-2)*fs error("Nyquist criterion is violated"); end  
    if fc * nfft <= fs error("fc is almost zero"); end

    % ---------------------------------------------------------------
    % Construction of frequency domain window function
    % ---------------------------------------------------------------
    
    winsize = 3;
    winsize = winsize - rem(winsize,2) + 1;
    
    mwind = fix(nfft/winsize);
    lby2  = (winsize - 1)/2;
    theta = -lby2:lby2;
    
    opwind = ones(winsize,1) * (theta .^2);
    opwind = opwind + opwind' + theta' * theta;
    opwind = 1 - (2*mwind/nfft)^2 * opwind;
    
    hex = ones(winsize,1) * theta;
    hex = abs(hex) + abs(hex') + abs(hex+hex');
    hex = (hex < winsize);
    
    opwind = opwind .* hex;
    opwind = opwind * (4 * mwind^2) / (7 * pi^2) ;

    % ---------------------------------------------------------------
    % Create a table to store bispectra and sleep stage labels
    % ---------------------------------------------------------------
    
    names = [X.Properties.VariableNames{channel}, "Annotations"];
    types = ["cell", "string"];
    
    bis = table(                        ...
        'Size',          [nrecs 2],     ...
        'VariableTypes', types,         ...
        'VariableNames', names);
    
    bis.Annotations = X.Annotations;

    % ---------------------------------------------------------------
    % Frequency axis
    % ---------------------------------------------------------------

    if (rem(nfft,2) == 0)
        freq = [-nfft/2:(nfft/2-1)];
    else
        freq = [-(nfft-1)/2:(nfft-1)/2];                        
    end

    % ---------------------------------------------------------------
    % Estimation of Bispectra
    % ---------------------------------------------------------------

    % array of indices for obtaining a truncated FFT
    idx = 1:nfft;                                                                                                            
    idx = idx(-fc * nfft <= freq * fs & freq * fs <= fc * nfft);
    len = numel(idx);

    % array of indices for accumulating triple products
    mask = hankel([1:len],[len,1:len-1] );

    % Bispectrum estimations for every 30sec record
    for i = 1:nrecs
        B  = zeros(len,len);                                              
        locseg = [1:nsamp]';
        x = cell2mat(X{i,channel});

        % Calculate partial bispectra for every segment
        for j = 1:K
            % Estimate the FFT of the segment
            xseg = x(locseg);                                              
            Xf   = fft(xseg-mean(xseg), nfft)/nsamp;
            
            % Discard any frequencies above fc and below -fc 
            % by truncating the result of the previous FFT
            Xf = fftshift(Xf);
            Xf = Xf(idx);                                                
            Xf = ifftshift(Xf);

            % Accumulate triple products
            CXf    = conj(Xf);
            B      = B + (Xf * Xf.') .* reshape(CXf(mask), len, len);    
            locseg = locseg + nadvance;                                    
        end

        % Shift the elements of the bispectrum matrix
        B = fftshift(B)/K;

        % Frequency domain smoothing
        B = conv2(B,opwind,'same');    
        
        % Save the results and proceed to the next record
        bis{i,1} = {B};
    end

    % ---------------------------------------------------------------
    % Truncate the frequency axis and rescale according to fs 
    % ---------------------------------------------------------------
    freq = freq(-fc * nfft <= freq *fs & freq * fs <= fc * nfft);
    freq = freq * fs / nfft;
return