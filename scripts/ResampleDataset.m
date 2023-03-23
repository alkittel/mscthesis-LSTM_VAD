clc
close all force
clear

%%
% load datastores.mat
% 
% outputLoc = '/Users/akittel/Documents/dutten/yr7 thesis/speech_corpus/trials_16k';
% 
% writeall(ads,outputLoc,'WriteFcn',@myWriter)
% 
% function myWriter(data,info,~)
%   fs = info.ReadInfo.SampleRate;
%   desiredFs = 16e3;
%   data = resample(data,desiredFs,fs);
%   data = data/(max(abs(data)));
%   audiowrite(info.SuggestedOutputName,data,desiredFs);
% end

%%
sAudioFolder='/Users/akittel/Documents/dutten/yr7 thesis/EasyComDataset/**/';
target='/Users/akittel/Documents/dutten/yr7 thesis/easycom_16k/';
eFiles=dir(fullfile(sAudioFolder,'*.wav'));
for ii=1:length(eFiles)
    sAudioFile=fullfile(eFiles(ii).folder,eFiles(ii).name);
    [y,Fs] = audioread(sAudioFile); y_resamp = resample(y,16e3,48000);
    sAudioFileOut=fullfile([target,extractAfter(eFiles(ii).folder,'EasyComDataset'),'/',eFiles(ii).name]);
    audiowrite(sAudioFileOut,y_resamp,16e3);
end

