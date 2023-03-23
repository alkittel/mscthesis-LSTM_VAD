close all force
clear
clc

T = 300;
parameters

load('perfeval_wsa.mat')
target = [labels{:}];

matFileObj = matfile('WSAdata.mat');
ads = matFileObj.test;
lss = matFileObj.lss;
reset(ads)


%% crosstalk shiz
tic
parfor ii = 1:numel(labels)
    disp(num2str(ii))
    bRef = sigroi2binmask(round( ...
        getLabelValues(lss,ii,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
        N);
    % crosstalk label
    bCt = false(size(bRef));
    crList = dir(['**/trials_16k/**/*',ads.Files{ii}(end-18:end-4),'*.wav']);
    mainIdx = str2double(extractBetween(string(ads.Files{ii}),'Talker','_'));    
    for hh = 1:numel(crList)
        lssIdx = find(ismember(lss.Source,[crList(hh).folder,'/',crList(hh).name]));
        crossIdx = str2double(extractBetween(crList(hh).name,'Talker','_'));    
        bTrial = sigroi2binmask(round( ...
            getLabelValues(lss,lssIdx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
            N);
        if crossIdx ~= mainIdx
            bCt = bCt | bTrial;
            disp('crosstalkers')
        end
    end
    crosstalk{ii,1} = logical(interp1(tSec,double(bCt&~bRef),tSecFrame,'nearest','extrap'));
end
toc

%% roc
disp(netNames.')
for ii = 1:size(NET_OUT,2)
    netCross = roc([crosstalk{:}],[NET_OUT{:,ii}]);
    netCrossHFA(ii,1) = netCross.TPR*100;
end
disp('Network crosstalk tpr:')
disp(netCrossHFA)

for ii = 1:size(VAD_OUT,2)
    vadCross = roc([crosstalk{:}],[VAD_OUT{:,ii}]);
    vadCrossHFA(ii,1) = vadCross.TPR*100;
end
disp('VAD crosstalk tpr:')
disp(vadCrossHFA)