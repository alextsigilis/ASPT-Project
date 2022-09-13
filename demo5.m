% ====================================================================
% Author: Christodoulos Michaelides
% Date: August, 2022
% --------------------------------------------------------------------
%
% Demo script 5:
%
% 1) Choose a number of patients from the dataset
% 2) Choose one out of the four EEG channels and decompose it to 
%    delta, theta, alpha and beta waves by applying the discrete 
%    wavelet transform.
% 3) Plot histograms of the DWT coefficients 
% ====================================================================

% ====================================================================
% 1) Reset your Workspace
% ====================================================================
clear all; close all; clc;

% ====================================================================
% 2) Script parameters (Choose any value you want to experiment with)
% ====================================================================

% n1: (integer 1-154) ID of first patient 
% n2: (integer 1-154) ID of last patient
% channel: (integer 1-4) ID of EEG channel 
n1 = 1; n2 = 25;
channel = 1;

% --------------------------------------------------------------------
% Do not change anything below that point, unless you know what you 
% you are doing.
% --------------------------------------------------------------------
names = {'delta' 'theta' 'alpha' 'beta' 'Annotations'};
types = {'cell'  'cell'  'cell'  'cell' 'string'};
X = table('Size',[0 5], 'VariableTypes', types, 'VariableNames',names);

% ===================================================================
% 3) Perform MRA decomposition in the selected EEG channel of every 
% patient
% ===================================================================

for i = n1:1:n2
    % Progress status
    fprintf("Processing patient %d\n\n",i);

    % Make sure the input file exists
    % before attempting to open it
    edf = sprintf("SN%03d.edf",i);
    mat = sprintf("%03d.mat",i);
    if (~isfile(edf)) && (~isfile(mat))
        fprintf("Patient %d does not exist\n\n",i);
        continue;
    end

    % Load the EEG recordings and sleep stage labels
    Z = loadEDF(i);

    % Estimate DWT coefficients for every frequency scale
    Y = mraEEG(Z,channel);

    % save the result and move on to the next patient
    X = vertcat(X, Y);
end

% ===================================================================
% 4) Plot histograms to visualize the distributions of the DWT
% coefficients in every sleep stage and MRA level.
% ===================================================================

% Array of valid EEG waves
waves = ["delta" "theta" "alpha" "beta"];

% Array of valid sleep stages
stages = [
    "Sleep stage W",        ...
    "Sleep stage R",        ...
    "Sleep stage N1",       ...
    "Sleep stage N2",       ...
    "Sleep stage N3"        ...
];

% mask: A colection of binary masks for different sleep stages
% mask{1}: binary mask for selecting "Sleep stage W" rows
% mask{2}: binary mask for selecting "Sleep stage R" rows
% mask{3}: binary mask for selecting "Sleep stage N1" rows
% mask{4}: binary mask for selecting "Sleep stage N2" rows
% mask{5}: binary mask for selecting "Sleep stage N3" rows
mask = cell(5,1);

for i = 1:1:numel(stages)
    mask{i} = (X.Annotations == stages(i));
end

for j = 1:1:numel(waves)
    figure(j); hold on;
    xlabel('DWT coefficients');
    ylabel('Relative frequency');
    title(sprintf('Histogram of %s wave coefficients', waves(j)));

    for i = 1:1:numel(stages)
        x = cell2mat(X{mask{i},j});
        [counts, centers] = hist(x, 1600);
        counts = counts / sum(counts);
        plot(centers, counts);
    end

    legend(stages);
end