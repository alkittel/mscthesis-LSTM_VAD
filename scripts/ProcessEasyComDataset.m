clear

cd('/Users/akittel/Documents/dutten/yr7 thesis')

addpath(genpath(pwd))

T = 60;
parameters


%% Generate recording datastores
%
%
% Recording and label paths
filePath = '/Users/akittel/Documents/dutten/yr7 thesis/easycom_16k/Close_Microphone_Audio';
labelPath = '/Users/akittel/Documents/dutten/yr7 thesis/easycom_16k/Speech_Transcriptions/';

% Generate datastores (audio datastore and labeled signal set
ads = audioDatastore(filePath,'FileExtensions','.wav', ...
    'IncludeSubfolders', true, 'LabelSource', 'foldernames');
lss = genEasycomLabels(ads,labelPath);
reset(ads)
kj = 0;
while hasdata(ads)
kj = kj+1;
    disp('another one')
    % Read recording and label files
    [sig,sigInfo] = read(ads);
    idx = find(ismember(lss.Source,sigInfo.FileName));
    blagh = sigroi2binmask(round( ...
        getLabelValues(lss,idx,'MainTalkerRegions').ROILimits * fsHz) + 1, ...
        numel(sig));
    ref{kj,1} = blagh.';
end

save('src/MAT-files/ECdataNEW.mat', ...
    "lss","ads","ref")

function lss = genEasycomLabels(ads,labelPath)
reset(ads)

% signal label definitions
dMainTalkerRegions = signalLabelDefinition('MainTalkerRegions', ...
    'LabelType','roi', ...
    'LabelDataType','logical', ...
    'Description','Regions of labeled main talker activity');

dTalkerID = signalLabelDefinition('TalkerID', ...
    'LabelType','attribute', ...
    'LabelDataType','numeric', ...
    'Description','Talker number ID');

lbldefs = [dMainTalkerRegions dTalkerID];

% create labeled signal set
lss = labeledSignalSet(ads,lbldefs,'SampleRate',16e3, ...
    'Description','Main talker activity ROIs');

% assign correct labels to corresponding recording depending on talker ID
kj = 1;
for ii = 1:672
    x = audioread(ads.Files{kj});

    % identify label file from recording time stamp
    timeStamp = ads.Files{kj}(end-29:end-21);
    jName = [labelPath,char(ads.Labels(kj)),'/',timeStamp,'.json'];

    % identify talker ID of current datastore member and session
    tID = str2double(ads.Files{kj}(end-4));

    % extract ROI times from label file corresponding to talker ID
    data = getECjsonData(tID,jName);
    starts  = [data.Start_Frame]/20;
    ends = [data.End_Frame]/20;

    mtr = [starts' ends'];
    mtrsz = [size(mtr,1) 1];

    % set talker LSS label values (speaker ROI and talker number)
    setLabelValue(lss,kj,'MainTalkerRegions', ...
        mtr,true(mtrsz));

    setLabelValue(lss,kj,'TalkerID', ...
        tID);

    bRef = sigroi2binmask(round( ...
        getLabelValues(lss,kj,'MainTalkerRegions').ROILimits * 16e3) + 1, ...
        numel(x));

    kj = kj + 1;
end

reset(ads)

end

function out = getECjsonData(tID,jName)
fid = fopen(jName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file

data = jsondecode(str); % Using the jsondecode function to parse JSON from string
out = data([data.Participant_ID] == tID); % Grab start and stop speak times for tID participant
end