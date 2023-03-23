function F = dynamicMeasures(target,test)
%dynamicMeasures    Calculate 4 dynamic VAD error types relative to a
%target reference (Freeman et al. 1988 [1]) from binary mask speaker
%activity
%
%USAGE
%   F = dynamicMeasures(target,test)
%
%INPUT ARGUMENTS
%   target : target utterance binary mask
%   test   : test utterance binary mask
%
%OUTPUT ARGUMENTS
%   F : Freeman et al. error structure
%       .FEC  - clipping at the front end of speech [ms]
%                   (time between reference onset of utterance and next
%                   test detection start)
%       .OVER - hangover after speech [ms]
%                   (time between target utterance offset and next test
%                   detection stop. set to NaN if the next stop falls
%                   within a target utterance ROI)
%       .MSC  - clipping in the middle of speech [%]
%                   (rate of 0's in utterances, disregarding FEC)
%       .NDS  - noise detected as speech [%]
%                   (rate of 1's outside of utterances, disregarding OVER)
%       .nC   - number of target utterances caught by at least one test
%                   detection
%
%REFERENCES
%   [1] - DK Freeman, CB Southcott, I Boyd, in Proc. of IEE
%         Colloquium on Digitized Speech Communication via Mobile Radio.
%         A voice activity detector for the Pan-European digital cellular
%         mobile telephone service (IEEE, London, United Kingdom, 1988)
%
%   Developed with Matlab 9.12.0.2039608 (R2022a)
%
%   Author :
%       Alexander Kittel (student, s164010)
%       Technical University of Denmark
%       Dept. of Health Technology
%
%   Date :
%       2022/11/07
%   ***********************************************************************


%% CHECK INPUTS
%
%
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check dimensionality
if size(target) ~= size(test)
    error('Target and test pattern must be of equal size.')
end

% Check if target and test are binary
if mean((target(:) == 0) + (target(:) == 1)) < 1 || ...
        mean((test(:) == 0) + (test(:) == 1)) < 1
    error('Target and test pattern must be binary.')
end


%% FREEMAN ET AL. DYNAMIC VAD METRIC ANALYSIS
%
%
% ROI tables
targetSec = binmask2sigroi(target)/16e3;
targetIdx = binmask2sigroi(target);
nUtterances = height(targetSec);

testSec = binmask2sigroi(test)/16e3;

% Prealloc
FEC  = zeros(nUtterances,1);
OVER = zeros(nUtterances,1);
MSC  = 0;
NDS  = 0;

% Init
backstop = 1;

% Loop over number of reference utterances
for ii = 1:nUtterances
    % FEC & OVER
    onset = targetSec(ii,1);
    offset = targetSec(ii,2);

    % Test detection starts and stops within utterance
    detectionStarts = ...
        testSec(testSec(:,1) >= onset & testSec(:,1) <= offset,1);

    detectionStops = ...
        testSec(testSec(:,2) >= onset & testSec(:,2) <= offset,2);
    if isempty(detectionStops)
        detectionStops(1) = nan;
    end

    % First immediate detection start following current utterance offset
    nextDetectionStart = testSec(find(testSec(:,1) > offset,1),1);

    % First immediate detection stop following current utterance offset
    nextDetectionStop = testSec(find(testSec(:,2) > offset,1),2);

    % FEC = positive distance from onset to next detection
    if isempty(detectionStarts) || (detectionStarts(1) > detectionStops(1))
        % if no detection starts in utterance, no FEC
        % if the first detection start is after a detection stop within
        % utterance, no FEC
        FEC(ii) = nan;
        nFEC = 0;
    else
        % if first detection start inside utterance
        FEC(ii) = detectionStarts(1) - onset;
        nFEC = round(FEC(ii)*16e3);
    end

    % OVER = positive distance from offset to next detection stop
    if target(end) == 1 && ii == nUtterances
        % overhang = 0 if signal ends within utterance
        OVER(ii) = 0;
        nOVER = round(OVER(ii)*16e3);
    elseif isempty(nextDetectionStop) || ii < nUtterances && targetSec(ii+1,1) < nextDetectionStop
        % case when detection spans 2 utterances or detections are ended
        OVER(ii) = nan;
        nOVER = 0;
    elseif nextDetectionStart > nextDetectionStop
        % overhang calculated only for detections within current utterance
        OVER(ii) = nextDetectionStop - offset;
        nOVER = round(OVER(ii)*16e3);
    else
        nOVER = 0;
    end

    % MSC & NDS
    nOnset = targetIdx(ii,1);
    nOffset = targetIdx(ii,2);

    % MSC = # of wrong detections (0) inside utterance disregarding FEC
    MSCrange = test(nOnset+nFEC:nOffset);
    MSC = MSC + sum(~MSCrange);

    % NDS = # of wrong detections (1) between backstop and this onset
    % disregarding OVER
    NDSrange = test(backstop:nOnset);
    NDS = NDS + sum(NDSrange);

    % new backstop = last offset + overhang
    backstop = nOffset + nOVER;
end


%% OUTPUT STRUCTURE
%
%
F.FEC  = FEC*1e3;
F.OVER = OVER*1e3;
F.MSC  = MSC/numel(target)*100;
F.NDS  = NDS/numel(target)*100;

end
