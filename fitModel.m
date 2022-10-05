% ===================================================================
% Author: Christodoulos Michaelides
% Date: October 5th, 2022
% -------------------------------------------------------------------
%
% Script Description:
%
% Using the extracted features as predictors, train a classifier
% for automatic Sleep-stage-scoring
% ===================================================================

clear all;
close all;
clc;

N = 20;
learner = "svm";

X_train = [];
y_train = [];

for i = 1:1:N
    filename = sprintf("%d.mat",i);
    
    if ~isfile(filename)
        continue;
    end

    load(filename);

    X_train = [X_train; X];
    y_train = [y_train; y];
end

model = fitcecoc(           ...
    X_train,y_train,        ...
    "Coding","onevsone",    ...
    "Learners",learner,     ... 
    "Verbose",2);