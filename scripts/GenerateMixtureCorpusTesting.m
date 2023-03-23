
clc
clear
rng(0)

%% paths
%
audioRoot = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/TIMIT/';
noiseRoot = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/DEMAND/test/';
audioList = dir([audioRoot,'TEST/**/*.wav']); % files from TIMIT train subset
noiseList = dir([noiseRoot,'**/*.wav']);

%% params
%
snrdB = [-10 -6 -3 0 2 4 6 10 15];
noiseTypes = ["kitchen","living","office"];
nNoiseTypes = numel(noiseTypes);
nSNRs = numel(snrdB);

nFiles = numel(audioList);
nMixtures = 50;

conditions  = cartesianProduct({1:nMixtures,1:nSNRs,1:nNoiseTypes});
nConditions = size(conditions,1);
idxMixtures = conditions(:,1);
idxSNR      = conditions(:,2);
idxNoise    = conditions(:,3);
snrdBVec    = snrdB(idxSNR);
noiseVec    = noiseTypes(idxNoise);

fsHz = 16e3;
T = 120;
N = T*fsHz;

%% mixtures
%
% gap - sentence - gap - sentence - gap - etc totalling N_WSA frames
% tic
disp('DJ KHALEEEEED')
parfor ii = 1:nConditions
    tic
    clean = zeros(N,1);
    bVAD = zeros(1,N);

    id = 0.5*fsHz;
    % break if within 20 seconds of recording duration
    while id < N - 20*fsHz
        % random TIMIT
        timit = audioList(randi(numel(audioList),1,1)); %#ok<*PFBNS> 
        x = audioread([timit.folder,'/',timit.name]);
        [~,~,vad,~] = labelsTIMIT([timit.folder,'/',timit.name]);

        % put TIMIT file at index
        clean(id:id+numel(x)-1) = x;
        bVAD(id:id+numel(vad)-1) = vad > -1;

        % jump ahead
        id = id + 5*numel(x);
    end

    %
    clean = clean(1:N);

    % mix in noise
    noise = audioread(noiseList(randi(numel(noiseList),1,1)).name);
    mix = adjustSNR(clean,noise(1:N),snrdBVec(ii));

    mixes{ii,1} = mix;
    labels{ii,1} = bVAD;

    disp('another one')
    toc
end
% toc
clear audioList audioRoot noiseList noiseRoot
save('/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/MIXtestdata.mat','-v7.3')