% Orfanos Dimitrios, 9579
%
% ==========================================================
% 
% Function Description: 
%
% This function decomposes an entire ECG recording into
% low, mid and high frequencies by applying the
% DWT algorithm.
% ----------------------------------------------------------
% Function Arguments:
%
% X: (timetable) the EEG recordings and sleep stage
% labels. You should use loadEDF() to create this
% timetable in order to ensure that its dimensions are
% compatible with this function.
%
% ----------------------------------------------------------
% Return Values:
%
% s: (table) A table with 1 + 1 columns.
% column 1: low frequencies DWT coefficients for every segment
% column 2: sleep stage annotations for every segment
%
% This is what the table should look like
%
%      low         Annotations     
% ______________  ______________ 
%
% {248×1 double} "Sleep stage W"  
% {248×1 double} "Sleep stage R"  
%       :               :                       
% {248×1 double} "Sleep stage W"  
% {248×1 double} "Sleep stage W"  
%
% ==========================================================
function s = mralevel2(X,channel,wavelet)
    % fs: (integer) sampling rate of X{:,idx} in Hertz
	% w:  (integer) duration of labeled segments in seconds
	% K:  (integer) number of labeled segments
	% N:  (integer) number of samples per labeled segment
	fs = 256; w = 30; N = w * fs; K = size(X,1); 
	
    % extract the EEG recording from the timetable
	x = X{:,channel};
	x = cell2mat(x);
	x = reshape(x, [N K]);

    % Initialize an empty table to store the DWT
    % coefficients and the sleep stage labels
    names = {'low' 'Annotations'}; 
    s  = cell(K,2);
    s  = cell2table(s,'VariableNames',names);

    % Decompose every 30sec segment into 3 levels
    % level 1 details:  64Hz - 128Hz    high frequencies
    % level 2 details:  32Hz - 64Hz     mid frequencies
    % level 3 details:  0Hz - 32Hz      low frequencies 

    for i = 1:K
        [A1,~] = dwt(x(:,i),wavelet);%A1 is 0-64
        [A2,~] = dwt(A1,wavelet);%A2 is 0-32
        A1 = idwt(A2,zeros(length(A2),1),wavelet);
        sig = idwt(A1,zeros(length(A1),1),wavelet);
        s.low{i} = sig(1:7680);
    end

    s.Annotations = X.Annotations;
end