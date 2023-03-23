clear

cd('/Users/akittel/Documents/dutten/yr7 thesis')

addpath(genpath(pwd))

tic


%% Generate recording datastores 
% 
%
% Recording and label paths
triadPath = 'speech_corpus/trials_16k/';
labelPath = 'speech_corpus/labels/';

% Generate datastores (audio datastore and labeled signal set
[ads,lss] = genCorpusDatastore(triadPath,labelPath);


%% Save datastores in mat file for future use
% 
% 
save('src/MAT-files/datastoresDS.mat', ...
    "lss","ads")

timeSec = toc;
disp(['Time elapsed: ' num2str(timeSec)])