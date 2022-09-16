% =====================================================================
% Author: Christodoulos Michaelides
% Date: September 7th, 2022
% ---------------------------------------------------------------------
%
% Function Description:
% This function can be used to detect quadratic phase coupling (QPC)
% between different types of EEG waves.
% ---------------------------------------------------------------------
%
% Arguments List: (b)
%
% b: (table) the bicoherence matrices of an EEG channel. You must 
% use bicEEG() to obtain this table.
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

function QPC = findQPC(b)   

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
    maxima = find(maxima == 1);

    % Convert linear indices to subscripts
    [x, y] = ind2sub(size(bic),maxima);
    qpc = bic(maxima); qpc = qpc(:);

    QPC{i,"peaks"} = {[x y qpc]};
end

% Copy sleep stage Annotations
QPC.Annotations = b.Annotations;

end