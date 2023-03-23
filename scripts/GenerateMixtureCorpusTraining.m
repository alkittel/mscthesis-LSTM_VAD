
clc
clear
rng(0)

%% paths
%
audioRoot = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/TIMIT/';
noiseRoot = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/DEMAND/train/';
audioList = dir([audioRoot,'TRAIN/**/*.wav']); % files from TIMIT train subset
noiseList = dir([noiseRoot,'PCAFETER/*.wav']);
target = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/mix-in training wavs/';


%% params
%
nNoiseTypes = 1;
snrdB = -20:5:20;
nSNRs = numel(snrdB);

nFiles = numel(audioList);
    nMixtures = 20;
    
    conditions  = cartesianProduct({1:nMixtures,1:nSNRs,1:nNoiseTypes});
    nConditions = size(conditions,1);
    snrdBVec    = snrdB(conditions(:,2));
    
    fsHz = 16e3;
    T = 300;
    N = T*fsHz;
    tSec = linspace(0,300,N);
    winType = 'hann';
    winSec = 15e-3;
    winSize = ceil(winSec * fsHz);
    overlapLength = winSize / 2;
    stepSize = winSize - overlapLength;
    nFrames = round((N-overlapLength)/stepSize);
    tSecFrame = (winSize/2 + (0:nFrames-1) * stepSize)/fsHz;


%% mixtures
%
% gap - sentence - gap - sentence - gap - etc totalling N_WSA frames
disp('DJ KHALEED')
parfor ii = 1:nConditions
    clean = zeros(N,1);
    bVAD = zeros(1,N);

    id = 0.5*fsHz;
    % break if within 20 seconds of recording duration
    while id < N - 20*fsHz
        % random TIMIT
        timit = audioList(randi(numel(audioList),1,1));
        x = audioread([timit.folder,'/',timit.name]);
        [~,~,vad,~] = labelsTIMIT([timit.folder,'/',timit.name]);

        % put TIMIT file at index
        clean(id:id+numel(x)-1) = x;
        bVAD(id:id+numel(vad)-1) = vad > -1;

        % jump ahead
        id = id + 5*numel(x);
    end

    % truncate
    clean = clean(1:N);

    labels{ii,1} = bVAD(1:N);

    % mix in cafeteria noise
    noise = audioread(noiseList(randi(numel(noiseList),1,1)).name);
    mix = adjustSNR(clean,noise(1:N),snrdBVec(ii));
    
    audiowrite([target,num2str(ii),'.wav'],rescale(mix,-1,1),fsHz)       

    disp('another one')
end

save([target,'labels.mat'], ...
    'labels','conditions','snrdBVec','-v7.3')