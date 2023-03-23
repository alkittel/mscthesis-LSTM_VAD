%% DOUBLE CHECK THAT THESE FEATURES ARENT ALREADY SAVED SOMEWHERE
%
%
clc
clear
rng(0)

T = 300; % WSA duration
parameters


%% PREPARE FEATURE EXTRACTION AND DATASET ORGANIZATION
%
%
% Load datastores of recordings and labels
load WSAMIXdata.mat train validate test lss

% Target path for saving feature dataset
featuresetsTargetPath = '/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/featuresets/';

% Feature dataset filename (NAME AFTER FEATURES EXTRACTED PLEASE)
featureFilename = 'WSAMIX_8MFs.mat';

% Feature time vector
tSecFeat = (winSize/2 + (0:nFrames-1) * stepSize)/fsHz;

% Recording time vector
tSecRef = linspace(0,T,N);

% Built-in feature extraction
afe = audioFeatureExtractor( ...
    SampleRate = fsHz, ...
    Window = feval(winType,winSize,'periodic'), ...
    OverlapLength = overlapLength, ...
    SpectralDescriptorInput = 'melSpectrum', ...
    melspectrum        = true);

% AFE parameters
setExtractorParameters(afe,'melSpectrum',NumBands=8)


%% EXTRACT FEATURES AND ALIGN LABELS
%
%
% Prealloc
XTrain = cell(numel(train.Files),1); % training features
YTrain = cell(numel(train.Files),1); % training labels

XValid = cell(numel(validate.Files),1); % validation features
YValid = cell(numel(validate.Files),1); % validation labels

XTest = cell(numel(test.Files),1); % test features
YTest = cell(numel(test.Files),1); % test labels

% Reset datastores to first entry
reset(train)
reset(validate)
reset(test)

% It's lights out and AWAY WE GO!
tic

% Generate training feature dataset
kj = 0;
while hasdata(train)
    kj = kj+1;

    % Read recording and label files
    [sig,sigInfo] = read(train);
    idx = find(ismember(lss.Source,sigInfo.FileName));
    bRef = sigroi2binmask(round( ...
        getLabelValues(lss,idx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
        numel(sig));

    % Extract features
    disp(['Extracting training features ', ...
        num2str(kj),' of ',num2str(numel(train.Files))])

    XTrain{kj} = extract(afe,sig)';

    % Temporally align binary reference mask to features
    YTrain{kj} = categorical(interp1(tSecRef,double(bRef),tSecFeat,'nearest','extrap'));

end

% Generate validation feature dataset
kj = 0;
while hasdata(validate)
    kj = kj+1;

    % Read recording and label files
    [sig,sigInfo] = read(validate);
    idx = find(ismember(lss.Source,sigInfo.FileName));
    bRef = sigroi2binmask(round( ...
        getLabelValues(lss,idx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
        numel(sig));

    % Extract features
    disp(['Extracting validation features ', ...
        num2str(kj),' of ',num2str(numel(validate.Files))])

    XValid{kj} = extract(afe,sig)';

    % Temporally align binary reference mask to features
    YValid{kj} = categorical(interp1(tSecRef,double(bRef),tSecFeat,'nearest','extrap'));

end

% Generate testing feature dataset
kj = 0;
while hasdata(test)
    kj = kj+1;

    % Read recording and label files
    [sig,sigInfo] = read(test);
    idx = find(ismember(lss.Source,sigInfo.FileName));
    bRef = sigroi2binmask(round( ...
        getLabelValues(lss,idx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
        numel(sig));

    % Extract features
    disp(['Extracting testing features ', ...
        num2str(kj),' of ',num2str(numel(test.Files))])

    XTest{kj} = extract(afe,sig)';

    % Temporally align binary reference mask to features
    YTest{kj} = categorical(interp1(tSecRef,double(bRef),tSecFeat,'nearest','extrap'));

end

timeFeatureExtraction = toc;


%% SAVE MAT-FILE
%
%
% Featuremap
fullFeatureMap = info(afe);

%
save([featuresetsTargetPath,featureFilename], ...
    'XTrain','YTrain','XTest','YTest','XValid','YValid', ...
    'fullFeatureMap','afe','-v7.3')

%
timeTotal = toc;
disp(['Elapsed time is ',num2str(timeTotal/60,3), ...
    ' minutes, for an average of ', num2str(timeFeatureExtraction/900,3), ...
    ' seconds per recording.'])