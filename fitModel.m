% ===================================================================
% Authors: Christodoulos Michaelides, Stergios Grigoriou
% Date: October 5th, 2022
% -------------------------------------------------------------------
%
% Script Description:
%
% Using the extracted features as predictors, train a classifier
% for automatic Sleep-stage-scoring
% Run this script inside the features folder.
% ===================================================================

clear
close all
clc
tic

learner = "svm";

X_train = [];
y_train = [];
stages =[ ...
    "Sleep stage W",     ...
    "Sleep stage N1",    ...
    "Sleep stage N2",    ...
    "Sleep stage N3",    ...
    "Sleep stage R"];
%% Load Features
for i = [1:13,15:35,37:63,65:97,99:134,136:154]
    filename = sprintf("%d.mat",i);
    
    if ~isfile(filename)
        continue;
    end

    load(filename);
    % Store 5-epoch features at each row t-2,t-1,..,t+2
    X_train = [X_train; X];
    y_train = [y_train; y];
end
load('ecg_features.mat')
X_train = [X_train,T{:,1:3}];
X_train(isnan(X_train)) = 0;
X_train(isinf(X_train)) = 0;
X_train = (X_train - mean(X_train))./std(X_train);
% PCA
% [c,~,~,explained,~] = pca(X_train);
% ind = find(cumsum(explained)>95,1);
% X_train = X_train*c(:,1:ind);
%% Train model
model = fitcecoc(           ...
    X_train,y_train,        ...
    "Coding","onevsone",    ...
    "Learners",learner,     ... 
    "Verbose",1,'ClassNames',stages,'Options',statset('UseParallel',true));
%% Validate model
CVmodel = crossval(model,'Kfold',5,'Options',statset('UseParallel',true));
%% Accuracy 
l = length(X_train);
y_hat = string(model.predict(X_train));
acc1 = sum(y_hat == y_train)/l;
fprintf('Accuracy on train set is %.3f\n',acc1)
acc2 = CVmodel.kfoldloss;
fprintf('5 fold cross validation accuracy is %.3f\n',acc2)
y_hat_CV = string(CVmodel.predict(X_train));
kappa = cohensKappa(y_train,y_hat_CV);
fprintf('5 fold cross validation kappa is %.3f\n',kappa)
toc
function kappa = cohensKappa(y, yhat)
    C = confusionmat(y, yhat); % compute confusion matrix
    n = sum(C(:)); % get total N
    C = C./n; % Convert confusion matrix counts to proportion of n
    r = sum(C,2); % row sum
    s = sum(C); % column sum
    expected = r*s; % expected proportion for random agree
    po = sum(diag(C)); % Observed proportion correct
    pe = sum(diag(expected)); % Proportion correct expected
    kappa = (po-pe)/(1-pe); % Cohen's kappa
end