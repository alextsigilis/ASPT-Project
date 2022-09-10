% =====================================================================
% Author: Christodoulos Michaelides
% Date: September 7th, 2022
% ---------------------------------------------------------------------
%
% Function Description:
% This function can be used to detect quadratic phase coupling (QPC)
% between different types of EEG waves.
%
% Since we are interested in 4 types of EEG waves (delta, theta,
% alpha and beta), we need to distinguish 9 different types of 
% non-linear interactions:
%
% 1) Interactions between delta and delta waves
% 2) Interactions between delta and theta waves
% 3) Interactions between delta and alpha waves
% 4) Interactions between delta and beta  waves
% 5) Interactions between theta and theta waves
% 6) Interactions between theta and alpha waves
% 7) Interactions between theta and beta  waves
% 8) Interactions between alpha and alpha waves
% 9) Interactions between alpha and beta  waves
%
% Our algorithm detects QPC by splitting the primary region of the 
% bicoherence matrix into 9 different sub-regions and locating 
% the position and value of the absolute maxima for those sub-regions
% In total, 9 frequency pairs are extracted from every bicoherence 
% matrix, one for every matrix partition.
% ---------------------------------------------------------------------
%
% Arguments List: (b, f)
%
% b: (table) the bicoherence matrices of an EEG channel. You must 
% use bicEEG() to obtain this table.
%
% f: (1D array) the frequency axes of the bicoherence matrices. Once 
% again, you should use bicEEG() to obtain this array.
% ---------------------------------------------------------------------
% 
% Return List: (QPC)
% 
% QPC: (table) A table with 9 + 1 columns. The first 9 columns
% correspond to the 9 bicoherence partitions. The 10th column holds
% the Sleep stage Annotations. Every cell of the first 9 columns 
% contains a triplet of the following form:
% {f1,f2,val}
% where:
%   -> (f1,f2) is the frequency pair in which QPC takes place.
%   -> val is b(f1,f2). The "strength" of the QPC phenomenon.
%
% This is what the QPC table should look like:
%
%  ____________________________________________________________________________________
% | deltaDelta | deltaTheta | deltaAlpha | .... | alphaAlpha | alphaBeta | Annotations |
% |____________|____________|____________|______|____________|___________|_____________|
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |  
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |      .     |     .      |     .      |      |     .      |     .     |      .      |
% |      .     |     .      |     .      |      |     .      |     .     |      .      |
% |      .     |     .      |     .      |      |     .      |     .     |      .      | 
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |    {1x3}   |  {[1x3]}   |  {[1x3]}   | .... |  {[1x3]}   |  {[1x3]}  |    string   |
% |____________|____________|____________|______|____________|___________|_____________|
%
% ---------------------------------------------------------------------

function QPC = findQPC(b,f)

% f0: DC frequency (in Hertz)
% f1: upper limit for delta waves (in Hertz)
% f2: upper limit for theta waves (in Hertz)
% f3: upper limit for alpha waves (in Hertz)
% f4: upper limit for beta  waves (in Hertz)
f0 =  0.0;
f1 =  4.0;
f2 =  8.0;
f3 = 16.0;
f4 = 32.0;

% delta: delta waves region on bicoherence matrix
% theta: theta waves region on bicoherence matrix
% alpha: alpha waves region on bicoherence matrix
% beta:  beta waves region on bicoherence matrix
delta = (f > f0) & (f <= f1);  
theta = (f > f1) & (f <= f2);
alpha = (f > f2) & (f <= f3);
beta  = (f > f3) & (f <= f4);

% deltaDelta: interactions between delta and delta waves
% deltaTheta: interactions between delta and theta waves
% deltaAlpha: interactions between delta and alpha waves
% deltaBeta:  interactions between delta and beta  waves
deltaDelta = delta & delta' & (f >= f') & (f + f' <= f4);
deltaTheta = theta & delta' & (f >= f') & (f + f' <= f4);    
deltaAlpha = alpha & delta' & (f >= f') & (f + f' <= f4);     
deltaBeta  = beta  & delta' & (f >= f') & (f + f' <= f4);      

% thetaTheta: interactions between theta and theta waves
% thetaAlpha: interactions between theta and alpha waves
% thetaBeta:  interactions between alpha and beta  waves
thetaTheta = theta & theta' & (f >= f') & (f + f' <= f4);
thetaAlpha = alpha & theta' & (f >= f') & (f + f' <= f4); 
thetaBeta  = beta  & theta' & (f >= f') & (f + f' <= f4); 

% alphaAlpha: interactions between alpha and alpha waves
% alphaBeta:  interactions between alpha and beta  waves
alphaAlpha = alpha & alpha' & (f >= f') & (f + f' <= f4);
alphaBeta  = beta  & alpha' & (f >= f') & (f + f' <= f4);

% N: number of 30sec epochs per EEG channel
% epsilon: a small positive constant to ensure numerical
% stability when performing floating point operations
N = size(b,1);
epsilon = 1e-5;

names = [           ...
    "deltaDelta",   ...
    "deltaTheta",   ...
    "deltaAlpha",   ...
    "deltaBeta",    ...
    "thetaTheta",   ...
    "thetaAlpha",   ...
    "thetaBeta",    ...
    "alphaAlpha",   ...
    "alphaBeta",    ...
    "Annotations"];

types = [           ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "cell",         ...
    "string"        ...
    ];

QPC = table(                        ...
    'Size',             [N 10],     ...
    'VariableNames',    names,      ...
    'VariableTypes',    types);

for i = 1:N
    % Extract the bicoherence array
    bic = cell2mat(b{i,1}) + epsilon;
    
    % QPC frequency pairs in delta-delta waves
    temp = bic .* deltaDelta;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"deltaDelta"} = {[x, y, val]};
    
    % QPC frequency pairs in delta-theta waves
    temp = bic .* deltaTheta;
    [val , qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"deltaTheta"} = {[x, y, val]};
    
    % QPC frequency pairs in delta-alpha waves
    temp = bic .* deltaAlpha;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"deltaAlpha"} = {[x, y, val]};
    
    % QPC frequency pairs in delta-beta waves
    temp = bic .* deltaBeta;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"deltaBeta"}  = {[x, y, val]};

    % QPC frequency pairs in theta-theta waves
    temp = bic .* thetaTheta;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"thetaTheta"} = {[x, y, val]};
    
    % QPC frequency pairs in theta-alpha waves
    temp = bic .* thetaAlpha;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"thetaAlpha"} = {[x, y, val]};
    
    % QPC frequency pairs in theta-beta waves
    temp = bic .* thetaBeta;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"thetaBeta"}  = {[x, y, val]};

    % QPC frequency pairs in alpha-alpha waves
    temp = bic .* alphaAlpha;
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"alphaAlpha"} = {[x, y, val]};
    
    % QPC frequency pairs in alpha-beta waves
    temp = bic .* alphaBeta;    
    [val, qpc] = max(temp(:));
    [x, y] = ind2sub(size(bic), qpc);
    QPC{i,"alphaBeta"}  = {[x, y, val]};
end

% Copy sleep stage Annotations
QPC.Annotations = b.Annotations;

end