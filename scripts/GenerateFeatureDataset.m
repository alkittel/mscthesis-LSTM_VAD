%% DOUBLE CHECK THAT THESE FEATURES ARENT ALREADY SAVED SOMEWHERE
%
%
clc
clear
rng(0)

% Watch these
T = 300;
parameters

%% PREPARE FEATURE EXTRACTION AND DATASET ORGANIZATION
%
%
% Load datastores of recordings and labels
load WSAMIXdata.mat train validate test lss

% Target path for saving feature dataset
featuresetsTargetPath = '/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/featuresets/';

% Feature dataset filename (NAME AFTER DATASET)
featureFilename = 'WSAMIX_fullfeatures.mat';

% Built-in feature extraction
afe = audioFeatureExtractor( ...
    SampleRate = fsHz, ...
    Window = feval(winType,winSize,'periodic'), ...
    OverlapLength = overlapLength, ...
    SpectralDescriptorInput = 'melSpectrum', ...
    spectralCentroid = true, ...
    spectralCrest = true, ...
    spectralEntropy =  true, ...
    spectralFlatness = true, ...
    spectralSpread =   true, ...
    harmonicRatio = true, ...
    zerocrossrate =    true, ...
    shortTimeEnergy = true);

% AFE parameters
setExtractorParameters(afe,'melSpectrum',NumBands=8)
setExtractorParameters(afe,'mfcc',NumCoeffs=8)


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

    ggFeats = graf_ghosh_features(sig,fsHz, ...
        'winSec',winSec,'overlap',overlapLength/winSize);
    afeFeats = extract(afe,sig)';

    XTrain{kj} = vertcat(afeFeats,ggFeats);

    % Temporally align binary reference mask to features
    YTrain{kj} = categorical(interp1(tSec,double(bRef),tSecFrame,'nearest','extrap'));

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

    ggFeats = graf_ghosh_features(sig,fsHz, ...
        'winSec',winSec,'overlap',overlapLength/winSize);
    afeFeats = extract(afe,sig)';

    XValid{kj} = vertcat(afeFeats,ggFeats);

    % Temporally align binary reference mask to features
    YValid{kj} = categorical(interp1(tSec,double(bRef),tSecFrame,'nearest','extrap'));

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

    ggFeats = graf_ghosh_features(sig,fsHz, ...
        'winSec',winSec,'overlap',overlapLength/winSize);
    afeFeats = extract(afe,sig)';

    XTest{kj} = vertcat(afeFeats,ggFeats);

    % Temporally align binary reference mask to features
    YTest{kj} = categorical(interp1(tSec,double(bRef),tSecFrame,'nearest','extrap'));

end

timeFeatureExtraction = toc;


%% SAVE MAT-FILE
%
%
% Featuremaps 
afeFeatureMap = info(afe);
[~,ggFeatureMap] = graf_ghosh_features(sig,fsHz, ...
    'winSec',winSec,'overlap',overlapLength/winSize);

fullFeatureMap = struct( ...
    'melSpectrum',[1 2 3 4 5 6 7 8], ...
    'mfcc',[9 10 11 12 13 14 15 16], ...
    'spectralEntropy', 17, ...
    'spectralSpread', 18, ...
    'zerocrossrate', 19, ...
    'longTermSignalVariability',20, ...
    'modulationPhaseDifference',21, ...
    'modulationPhaseDifferenceHold',22, ...
    'combinedGrafModulationFeature',23, ...
    'shortTimeLogEnergy',24);

%
save([featuresetsTargetPath,featureFilename], ...
    'XTrain','YTrain','XTest','YTest','XValid','YValid', ...
    'fullFeatureMap','test','afe','-v7.3')

%
timeTotal = toc;
disp(['Elapsed time is ',num2str(timeTotal/60,3), ...
    ' minutes, for an average of ', num2str(timeFeatureExtraction/900,3), ...
    ' seconds per recording.'])