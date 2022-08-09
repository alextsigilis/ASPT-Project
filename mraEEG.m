% ==========================================================
% Author: Christodoulos Michaelides
% Date: August 9th, 2022
% ----------------------------------------------------------
% Function Description: 
%
% This function decomposes an entire EEG recording into
% delta, theta, alpha and beta waves by applying the
% DWT algorithm.
% ----------------------------------------------------------
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
% Return Values:
%
% coef: (table) A table with 4 + 1 columns.
% column 1: delta wave DWT coefficients for every segment
% column 2: theta wave DWT coefficients for every segment
% column 3: alpha wave DWT coefficients for every segment
% column 4: beta wave DWT coefficients for every segment
% column 5: sleep stage annotations for every segment
%
% This is what the table should look like
%
%     delta            theta           alpha           beta         Annotations  
% ______________  ______________  ______________  ______________  _______________
%
% {248×1 double}  {248×1 double}  {488×1 double}  {967×1 double}  "Sleep stage W"
% {248×1 double}  {248×1 double}  {488×1 double}  {967×1 double}  "Sleep stage R"
%       :               :               :               :                :       
% {248×1 double}  {248×1 double}  {488×1 double}  {967×1 double}  "Sleep stage W"
% {248×1 double}  {248×1 double}  {488×1 double}  {967×1 double}  "Sleep stage W"
%
% ==========================================================
function coef = mraEEG(X,idx)
    % fs: (integer) sampling rate of X{:,idx} in Hertz
	% w:  (integer) duration of labeled segments in seconds
	% K:  (integer) number of labeled segments
	% N:  (integer) number of samples per labeled segment
	fs = 256; w = 30; N = w * fs; K = size(X,1); 
	
    % extract the EEG recording from the timetable
	x = X{:,idx};
	x = cell2mat(x);
	x = reshape(x, [N K]);

    % Initialize an empty table to store the DWT
    % coefficients and the sleep stage labels
    names = {'delta' 'theta' 'alpha' 'beta' 'Annotations'}; 
    coef  = cell(K,5);
    coef  = cell2table(coef,'VariableNames',names);

    % Decompose every 30sec segment into 5 levels
    % level 1 details:  64Hz - 128Hz
    % level 2 details:  32Hz - 64Hz   gamma waves
    % level 3 details:  16Hz - 32Hz   beta waves
    % level 4 details:   8Hz - 16Hz   alpha waves 
    % level 5 details:   4Hz - 8Hz    theta waves
    % level 5 approx:    0Hz - 4Hz    delta waves
    for i = 1:1:K
        [C, L] = wavedec(x(:,i),5,'db5');
        coef.delta{i} = appcoef(C,L,'db5',5);
        [coef.theta{i},coef.alpha{i},coef.beta{i}] = detcoef(C,L,[5 4 3]);
    end

    coef.Annotations = X.Annotations;
end