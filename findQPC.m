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
% Our algorithm detects QPC with a simple 2-step process:
% 1) The local maxima of the bicoherence are located by using the
%    builtin imregionalmax() function of the Image Processing 
%    Toolbox directly on the bicoherence matrix.
% 
% 2) Local Maxima which do not exceed a user-specified threshold 
%    value are discarded.
%
% The remaining local maxima of the last step are considered 
% statistically significant.
% ---------------------------------------------------------------------
%
% Arguments List: (b, f)
%
% b: (table) the bicoherence matrices of an EEG channel. You must 
% use bicEEG() to obtain this table.
%
% f: (1D array) the frequency axes of the bicoherence matrices. Once 
% again, you should use bicEEG() to obtain this array.
%
% epsilon: (float in [0 1]) a threshold value to get rid of 
% insignificant local maxima
% ---------------------------------------------------------------------
% 
% Return List: (QPC)
% 
% QPC: (table) A table with 2 columns. 
% The first column contains a list (cell array) of triplets. Every 
% triplet contains the coordinates (f1,f2) of a local maximum 
% expressed in Hertz and the bicoherence value of that local maximum. 
% Therefore, the triplets will always have the following format:
% {[f1, f2, b(f1,f2)])
% The second column contains the Sleep stage Annotation for every 
% 30sec epoch.
% ---------------------------------------------------------------------

function QPC = findQPC(b,f,epsilon)   

% N: (integer) number of 30sec epochs per EEG channel
N = size(b,1);

% QPC: (table) This is where the bicoherence peaks are stored
QPC = table(                                    ...
    'Size',          [N 2],                     ...
    'VariableNames', ["peaks" "Annotations"],   ...
    'VariableTypes', ["cell" "string"]);

for i = 1:N
    % Extract the bicoherence array
    bic = cell2mat(b{i,1});

    % Locate all local maxima
    maxima = imregionalmax(bic);

    % Discard local maxima based on the
    % threshold imposed by epsilon
    maxima = maxima & bic >= epsilon;
    maxima = find(maxima == 1);

    % Convert linear indices to subscripts
    [x, y] = ind2sub(size(bic),maxima);
    qpc = bic(maxima); qpc = qpc(:);

    QPC{i,"peaks"} = {[x y qpc]};
end

% Copy sleep stage Annotations
QPC.Annotations = b.Annotations;

end