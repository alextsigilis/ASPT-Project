% ==========================================================
% Author: Christodoulos Michaelides
% Date: August 5th, 2022
% ----------------------------------------------------------
%
% Function Description: 
%
% This function decomposes an entire EEG recording into
% delta, theta, alpha and beta waves. It then estimates
% the variance, skewness and kurtosis of every wave type
% in a 30sec sliding window. Those estimates can then 
% be used by a classifier to identify the sleep stage 
% at every 30sec segment of the EEG.
% ----------------------------------------------------------
% 
% Function Arguments:
%
% X: (timetable) the EEG recordings and sleep stage
% labels. You should use loadEDF() to create this
% timetable in order to ensure that its dimensions are
% compatible with this function.
%
% idx: (integer) an index between 1 and 4. This index
% is used to select an EEG channel to decompose
% according to the following rule:
% idx = 1 => EEG F4-M1 is selected for analysis
% idx = 2 => EEG C4-M1 is selected for analysis
% idx = 3 => EEG O2-M1 is selected for analysis
% idx = 4 => EEG C3-M2 is selected for analysis
% ----------------------------------------------------------
%
% Return Values: delta, theta, alpha, beta
%
% delta: (table) a table with 3 + 1 columns. 
% The 1st column contains the standard deviation of
% the delta waves for each 30sec segment.
% The 2nd column contains the skewness of the delta
% waves for each 30sec segment.
% The 3rd column contains the kurtosis of the delta
% waves for each 30sec segment.
% The 4th column contains the sleep stage annotation
% for each 30sec segment
%
% theta: (table) same as delta but for theta waves instead
% 
% alpha: (table) same as delta but for alpha waves instead
%
% beta: (table) same as delta but for beta waves instead
% ----------------------------------------------------------
%
% Example: 
% Read the EEG recordings of patient No. 42, choose the 2nd 
% EEG channel, perform MRA decomposition and estimate higher  
% order statistics for each frequency scale.
%
% >> Z = loadEDF(42);
% >> [delta, theta, alpha, beta] = mraEEG(Z,2);
% >> disp(delta);
%
%     var         skw        krt        Annotations  
%    ______    _________    ______    _______________
%
%    178.59       1.1808    7.8672    "Sleep stage W"
%    62.976      0.72462    5.6213    "Sleep stage W"
%      .            .          .             .
%      .            .          .             .
%    8.3264     -0.27104    3.0522    "Sleep stage W"
% ==========================================================

function [delta, theta, alpha, beta] = mraEEG(X,idx)
    % fs: (integer) sampling rate of X{:,idx} in Hertz
	% w:  (integer) duration of labeled segments in seconds
	% K:  (integer) number of labeled segments
	% N:  (integer) number of samples per labeled segment
	fs = 256; w = 30; N = w * fs; [K,~] = size(X); 
	 
	x = X{:,idx};
	x = cell2mat(x);
	x = reshape(x, [N K]);		
	
    % wt:  array of wavelet coefficients
    % mra: array of wavelet frequency scales
	wt  = modwt(x,"db2",5);
	mra = modwtmra(wt,"db2");
		
    % delta: array of statistics for delta waves
    % theta: array of statistics for theta waves
    % alpha: array of statistics for alpha waves
    % beta:  array of statistics for beta waves
	delta = nan(K,3);
    theta = nan(K,3);
	alpha = nan(K,3);
	beta  = nan(K,3);
	
    % Estimate delta wave statistics
	delta(:,1) = std(mra(6,:,:),0,2);
	delta(:,2) = skewness(mra(6,:,:),0,2);
	delta(:,3) = kurtosis(mra(6,:,:),0,2);
	
    % Estimate theta wave statistics
	theta(:,1) = std(mra(5,:,:),0,2);
	theta(:,2) = skewness(mra(5,:,:),0,2);
	theta(:,3) = kurtosis(mra(5,:,:),0,2);
	
    % Estimate alpha wave statistics
	alpha(:,1)  = std(mra(4,:,:),0,2);
	alpha(:,2)  = skewness(mra(4,:,:),0,2);
	alpha(:,3)  = kurtosis(mra(4,:,:),0,2);
	
    % Estiamte beta wave statistics
	beta(:,1) = std(mra(3,:,:),0,2);
	beta(:,2) = skewness(mra(3,:,:),0,2);
	beta(:,3) = kurtosis(mra(3,:,:),0,2);
	
    % Variable names for each column
	names = {'var' 'skw' 'krt'};
	
    % Convert arrays into tables and rename their columns
	delta = array2table(delta, 'VariableNames', names);
	theta = array2table(theta, 'VariableNames', names);
	alpha = array2table(alpha, 'VariableNames', names);
	beta  = array2table(beta,  'VariableNames', names);
	
    % Add an extra column for sleep stage annotations
	delta = addvars(delta, X.Annotations, 'NewVariableNames', 'Annotations');
	alpha = addvars(alpha, X.Annotations, 'NewVariableNames', 'Annotations');
	beta  = addvars(beta,  X.Annotations, 'NewVariableNames', 'Annotations');
	theta = addvars(theta, X.Annotations, 'NewVariableNames', 'Annotations');
end