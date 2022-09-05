% =============================================================
% Author: Christodoulos Michaelides
% Date: August 22nd, 2022
% -------------------------------------------------------------
%
% Description:
% Reading EDF/EDF+ files seems to be extremely inefficient
% when using the builtin edfread() and edfinfo() functions.
% In our case things are made worse because of the fact that 
% patients have two EDF files associated with them. The first
% contains the EEG recordings and the second contains the sleep 
% stage Annotations. 
% 
% On the other hand, loading data from .mat files is 
% significantly faster. This script reads EEG recordings and 
% sleep stage Annotations from the EDF files and stores them
% in .mat files to speed up file I/O operations in the future.
%
% You don't have to manually load data from those .mat files. 
% Instead, loadEDF() attempts to retrieve data from them before
% reading the slower EDF files.
% =============================================================

clear all; close all; clc;

start = 1;            % first patient
stop  = 154;          % last patient

for i = start:1:stop
    matFile = sprintf("%03d.mat",i);
    edfFile = sprintf("SN%03d.edf",i);

    if ~isfile(edfFile)
        continue;
    end

    % Delete pre-existing mat-files
    if isfile(matFile)
        delete(matFile);
    end

    % Create new mat-files
    fprintf("Patient %d out of %d\n", i, stop);
    Z = loadEDF(i);
    save(matFile, 'Z');
end

clear all;
