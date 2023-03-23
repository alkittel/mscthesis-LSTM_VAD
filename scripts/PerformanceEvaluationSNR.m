%% SNR IS SO COOL AND ARBITRARY
%
%
%#ok<*PFBNS>

close all force
clear
clc

T = 120;
parameters


%% init

% mixes
mixFileObj = matfile("MIXtestdata.mat");
conditions = mixFileObj.conditions;
snrdB = mixFileObj.snrdB;

nConditions = size(conditions,1);

% Trained classifier
dirlist = dir('**/test/*.mat');
parfor ii = 1:numel(dirlist)
    disp(['Parpool working: ',num2str(ii),' of ',num2str(numel(dirlist))])
    mat = matfile(dirlist(ii).name);
    nets{ii} = mat.net;
    mu{ii} = mat.mu; sd{ii} = mat.sd;
    netNames{ii} = dirlist(ii).name(1:end-24);
    nMFCCs{ii} = str2double(strrep(netNames{ii}(end-16:end-15),'_',''));
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


%% classify
tic
parfor ii = 1:nConditions
    disp(['Working. ',num2str(ii),' out of ',num2str(nConditions),'.'])
    mix = mixFileObj.mixes(ii,1);
    ref = mixFileObj.labels(ii,1);
    x = mix{:};
    labels{ii} = logical(interp1(tSec,double(ref{:}),tSecFrame,'nearest','extrap'));
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
    end

    for kk = 1:nVads
        strMethod = ['vad',upper(methods{kk}{1}(1)),...
            lower(methods{kk}{1}(2:end))];
        % VAD
        [test,tSecVad] = feval(strMethod,x,fsHz,'initSec',0.5,...
            methods{kk}{2}{:});
        test = logical(interp1(tSecVad,double(test),tSecFrame,'nearest','extrap'));
        VAD_OUT{ii,kk} = test;
    end
end

toc

%% save
%
save('/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/results/final/perfeval_snr.mat', ...
    'NET_OUT','NET_SCR','VAD_OUT','labels','conditions','snrdB','netNames','dirlist','-v7.3')