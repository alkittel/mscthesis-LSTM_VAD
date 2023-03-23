clc
clear
close all


%%
load('/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/mix-in training wavs/labels.mat');
load('/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/WSAdata.mat')
adsWSA = ads;
lssWSA = lss;
clear ads lss

%%
recPath = '/Users/akittel/Documents/dutten/yr7 thesis/mixture_corpus/mix-in training wavs/';

adsMIX = audioDatastore(recPath,'FileExtensions','.wav', ...
    'IncludeSubfolders', true, 'LabelSource', 'foldernames');

dMainTalkerRegions = signalLabelDefinition('MainTalkerRegions', ...
    'LabelType','roi', ...
    'LabelDataType','logical', ...
    'Description','Regions of labeled main talker activity');

dTriadID = signalLabelDefinition('TriadID', ...
    'LabelType','attribute', ...
    'LabelDataType','string', ...
    'Description','Triad number ID');

lbldefs = [dMainTalkerRegions dTriadID];

% create labeled signal set
lssMIX = labeledSignalSet(adsMIX,lbldefs,'SampleRate',16e3, ...
    'Description','Main talker activity ROI');

reset(adsMIX)

while hasdata(adsMIX)
    [~,stuff] = read(adsMIX);
    midx = str2double(stuff.FileName(end-6:end-4));

    mtr = binmask2sigroi(logical(labels{midx}))/16e3;
    mtrsz = [size(mtr,1) 1];

    setLabelValue(lssMIX,midx,'MainTalkerRegions',mtr,true(mtrsz));
    setLabelValue(lssMIX,midx,'TriadID','synthetic');
    disp('we love fortnite')
end

[keep,toss] = splitEachLabel(adsWSA,.8,'randomized');
ads = audioDatastore(cat(1,keep.Files,adsMIX.Files),'LabelSource', 'foldernames');

lss = merge(lssWSA,lssMIX);


save('src/MAT-files/WSAMIXdata.mat', ...
    "lss","ads")