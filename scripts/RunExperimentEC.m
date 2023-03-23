%% INITIALIZATION
%
%
%#ok<*PFBNS>
%#function SeriesNetwork

clc
clear

cd('/Users/akittel/Documents/dutten/yr7 thesis')
addpath(genpath(pwd))

T = 60; % EC
parameters

rng(0)

%% Easycom (META) set
matFileObj = matfile('ECdata.mat');
ads = matFileObj.ads;
lss = matFileObj.lss;
reset(ads)

%% init

% Trained classifier
dirlist = dir('**/test/*.mat');
parfor ii = 1:numel(dirlist)
    disp(['Parpool working: ',num2str(ii),' of ',num2str(numel(dirlist))])
    mat = matfile(dirlist(ii).name);
    nets{ii} = mat.net;
    mu{ii} = mat.mu; sd{ii} = mat.sd;
    netNames{ii} = dirlist(ii).name(1:end-24);
end

nNets = numel(nets);

afeCell{1} = audioFeatureExtractor( ...
        SampleRate = fsHz, ...
        Window = feval(winType,winSize,'periodic'), ...
        OverlapLength = overlapLength, ...
        SpectralDescriptorInput = 'melSpectrum', ...
        mfcc=true);
afeCell{2} = audioFeatureExtractor( ...
        SampleRate = fsHz, ...
        Window = feval(winType,winSize,'periodic'), ...
        OverlapLength = overlapLength, ...
        SpectralDescriptorInput = 'melSpectrum', ...
        melSpectrum=true);
afeCell{3} = audioFeatureExtractor( ...
        SampleRate = fsHz, ...
        Window = feval(winType,winSize,'periodic'), ...
        OverlapLength = overlapLength, ...
        SpectralDescriptorInput = 'melSpectrum', ...
        mfcc=true);

setExtractorParameters(afeCell{1,1},'mfcc',NumCoeffs=8);
setExtractorParameters(afeCell{1,2},'melSpectrum',NumBands=8);
setExtractorParameters(afeCell{1,3},'mfcc',NumCoeffs=8);


% Other classifiers (Ghosh & Graf)
methods = {
    {'ghosh'  {'alpha' 0.6 'winSec' winSec}  }   ...
    {'graf'    {'alpha' 0.6 'winSec' winSec}  }   ...
    };
nVads = numel(methods);


%% RUN EXPERIMENT
%
%
tic
parfor ii = 1:numel(ads.Files)
    disp(['Working. ',num2str(ii),' out of ',num2str(numel(ads.Files)),'.'])
    x = audioread(ads.Files{ii});

    % run trained classifiers
    for jj = 1:nNets
        % features
        EXTRACTOR = afeCell{jj};
        feats = extract(EXTRACTOR,x)';
        % standardization
        feats2 = (feats-mu{jj})./sd{jj};
        % classification
        [test,score] = classify(nets{jj},feats2);
        NET_OUT{ii,jj} = test == '1';
        NET_SCR{ii,jj} = score(2,:); % true label score
        % label
        idx = find(ismember(lss.Source,ads.Files{ii}));
        bRef = sigroi2binmask(round( ...
            getLabelValues(lss,idx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
            numel(x));
        labels{ii} = logical(interp1(tSec,double(bRef),tSecFrame,'nearest','extrap'));
    end

    for kk = 1:nVads
        strMethod = ['vad',upper(methods{kk}{1}(1)),...
            lower(methods{kk}{1}(2:end))]
        % VAD
        [test,tSec_DUMB] = feval(strMethod,x,fsHz,'initSec',0,...
            methods{kk}{2}{:});
        test = logical(interp1(tSec_DUMB,double(test),tSecFrame,'nearest','extrap'));
        VAD_OUT{ii,kk} = test;
    end
    disp(['Done. ',num2str(ii),' out of ',num2str(numel(ads.Files)),'.'])
end
toc

%%
save('/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/results/final/perfeval_ec.mat', ...
    'NET_OUT','NET_SCR','VAD_OUT','labels','netNames','dirlist','-v7.3')