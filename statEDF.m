% =========================================================
% Author: Christodoulos Michaelides
% Date: August 1st, 2022
% ---------------------------------------------------------
%
% Function Description:
% This function estimates stanard deviation, skewness and
% kurtosis in a 30sec long time window for every available
% EEG channel in a timetable. (in our case there are 4
% EEG channels per patient)
% ---------------------------------------------------------
% 
% Function Arguments:
% Z: (timetable) a timetable which contains the recordings 
% of a single patient. Use loadEDF() to get the timetable.
% ---------------------------------------------------------
% 
% Return Variables:
% var, skw, krt: (tables). Every table consists of four
% columns (1 column per EEG channel) plus an extra column 
% for sleep stage annotations. 
% Every entry on these four columns contains the estimate
% of standard deviation, skewness or kurtosis for its 
% respective 30sec labeled segment.
%
% var:
%
% EEGF4_M1 | EEGC4_M1 | EEGO2_M1 | EEGC3_M2 | Annotations
% -------------------------------------------------------
%   var1   |   var2   |   var3   |   var4   |   label1  |
%   var5   |   var6   |   var7   |   var8   |   label2  |
%    .     |    .     |    .     |    .     |     .     |
%    .     |    .     |    .     |    .     |     .     |
% -------------------------------------------------------
%
% skw: 
%
% EEGF4_M1 | EEGC4_M1 | EEGO2_M1 | EEGC3_M2 | Annotations
% -------------------------------------------------------
%   skw1   |   skw2   |   skw3   |   skw4   |   label1  |
%   skw5   |   skw6   |   skw7   |   skw8   |   label2  |
%    .     |    .     |    .     |    .     |     .     |
%    .     |    .     |    .     |    .     |     .     |
% -------------------------------------------------------
%
% krt: 
%
% EEGF4_M1 | EEGC4_M1 | EEGO2_M1 | EEGC3_M2 | Annotations
% -------------------------------------------------------
%   krt1   |   krt2   |   krt3   |   krt4   |   label1  |
%   krt5   |   krt6   |   krt7   |   krt8   |   label2  |
%    .     |    .     |    .     |    .     |     .     |
%    .     |    .     |    .     |    .     |     .     |
% -------------------------------------------------------

function [var, skw, krt] = statEDF(Z)
    % M: (integer) number of EEG channels per patient
    % K: (integer) number of sleep stage labels
    [K,~] = size(Z); M = 4;
    var = nan(K,M);
    skw = nan(K,M);
    krt = nan(K,M);

    % go through every 30sec epoch
    for i = 1:1:K
        z = cell2mat(Z{i,1:M});
        var(i,:) = std(z,0,1);
        skw(i,:) = skewness(z,0,1);
        krt(i,:) = kurtosis(z,0,1);
    end

    % convert var,skw and
    % krt arrays into tables
	names = Z.Properties.VariableNames;
    var = array2table(var,'VariableNames',{names{1:M}});
	skw = array2table(skw,'VariableNames',{names{1:M}});
	krt = array2table(krt,'VariableNames',{names{1:M}});

    % add an extra column with sleep
    % stage annotations for every epoch
    var = addvars(var, Z.Annotations, 'NewVariableNames', 'Annotations');
    skw = addvars(skw, Z.Annotations, 'NewVariableNames', 'Annotations');
    krt = addvars(krt, Z.Annotations, 'NewVariableNames', 'Annotations');
end