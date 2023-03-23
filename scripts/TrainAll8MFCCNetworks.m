%% INITIALIZATION
%
%
clc
clear

cd('/Users/akittel/Documents/dutten/yr7 thesis')
addpath(genpath(pwd))

T = 300;
parameters

rng(0)

% Load feature dataset SMARTLY
disp('Loading features.mat ...')
tic
matObj = matfile('WSAMIX_8MFCCs.mat');
XTrain = matObj.XTrain;
YTrain = matObj.YTrain;
disp('Data loaded.')
toc


%% SELECT DESIRED FEATURES
%
%
% Standardize
mu = mean([XTrain{:}],2);
sd = std([XTrain{:}],[],2);

XTrain = cellfun(@(x)(x-mu)./sd,XTrain,'UniformOutput',false);

disp('Data standardized with means and std. deviations from training subset.')

%%
%% RESHAPE DATA
%
%
sequenceLength = 200;
sequenceOverlap = 0;

XTrain2 = ...
    helperFeatureVector2Sequence([XTrain{:}],sequenceLength,sequenceOverlap);

YTrain2 = ...
    helperFeatureVector2Sequence([YTrain{:}],sequenceLength,sequenceOverlap);


%% LSTM NETWORK ARCHITECTURE
%
%
layers = [ ...
    sequenceInputLayer(height(XTrain2{1}))

    batchNormalizationLayer

    lstmLayer(100,'OutputMode','sequence')

    dropoutLayer(0.2)

    lstmLayer(100,'OutputMode','sequence')

    dropoutLayer(0.2)

    fullyConnectedLayer(2)
    softmaxLayer
    classificationLayer];

disp(layers)

options = trainingOptions('adam');

options.MaxEpochs = 10;
options.MiniBatchSize = 32;
options.Shuffle = 'every-epoch';

options.Plots = 'none';
options.Verbose = true;


%% TRAIN NETWORK
%
%
[net,netInfo] = trainNetwork(XTrain2,YTrain2,layers,options);


%% SAVE DATA
%
%
% Save network info, standardization statistics and feature map
save(['src/MAT-files/trained VADs/test/2LSTM_BATCHNORM_DROPOUT20_8MFCCs_200SEQ_WSAMIX', ...
    datestr(now,'dd mmm yyyy HH MM SS'),'.mat'], ...
    'net','netInfo','options','mu','sd')
