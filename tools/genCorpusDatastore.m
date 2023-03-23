function [ads,lss] = genCorpusDatastore(triadPath,labelPath)
%genCorpusDatastore
%
% INPUTS:
%       triadPath - path to triad recordings 
%       labelPath - path to talker regions label files (.txt)
%      
% OUTPUTS:
%       ads - audio datastore, labeled by recording condition
%       lss - labeled signal set:
%           MainTalkerRegions - ROI of main talker regions in samples
%           TriadNumber - Triad number (1 to 25) of recording
%
% SEE ALSO:
%       WSA triad recording protocol


% sample rate
fsHz = 44.1E3;

% load audio datastore
ads = audioDatastore(triadPath,'FileExtensions','.wav', ...
    'IncludeSubfolders', true, 'LabelSource', 'foldernames');

% file label definition (recording condition, quiet as default)
ads.Labels(:) = categorical({'quiet'});

% signal label definitions
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
lss = labeledSignalSet(ads,lbldefs,'SampleRate',fsHz, ...
    'Description','Main talker activity ROIs and triad ID');

% assign correct labels to corresponding recording depending on talker ID
kj = 1;
while hasdata(ads)
    [~,info] = read(ads);

    % identify label file from recording time stamp
    timeStamp = info.FileName(end-11:end-4);
    txtName = dir([labelPath,'/**/*',timeStamp,'.txt']).name;

    % identify talker ID of current datastore member
    talker = str2double(extractBetween(info.FileName,'Talker','_Trial'));

    % identify triad number
    tn = extractBefore(txtName,'_');

    % extract ROI times from label file corresponding to talker ID
    mtr = genRefTalkerROI(txtName,talker);
    mtrsz = [size(mtr,1) 1];

    % set talker LSS label values (speaker ROI and triad number)
    setLabelValue(lss,kj,'MainTalkerRegions', ...
        mtr,true(mtrsz));
    setLabelValue(lss,kj,'TriadID', ...
        tn);

    % set recording condition (omni & dir are noisy, otherwise quiet) and
    % triad number as ADS labels
    if count(txtName,'Omni') || count(txtName,'Dir')
        ads.Labels(kj) = categorical({'noisy'});
    else
    end

    kj = kj + 1;
end

reset(ads)

end