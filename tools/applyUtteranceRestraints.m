function out = applyUtteranceRestraints(bVAD,uttSec,gapSec,stepSize,fsHz)
%applyUtteranceRestraints Applies minimum utterance and gap duration
%restrictions to a binary VAD output
%
%   Developed with Matlab 9.12.0.2039608 (R2022a)
%
%   Author :
%       Alexander Kittel (student, s164010)
%       Technical University of Denmark
%       Dept. of Health Technology
%
%   Date :
%       2022/11/13
%
%   Version  : v0.1.0    2022/11/13
%              v0.1.1    2022/11/18
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 5 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check if test is logical array
if ~islogical(bVAD)
    error('Target and test pattern must be binary.')
end

out = bVAD;


%% PREPROCESS
%
%
% Convert to ROI limits
roi = binmask2sigroi(bVAD);


%% APPLY MINIMUM DURATIONS FOR PAUSES
%
%
% Extract gap lengths in frames
gapFrames = zeros(1,height(roi));
gapFrames(1) = (roi(2,1) - roi(1,2));
for ii = 2:height(roi)
    gapFrames(ii) = (roi(ii,1) - roi(ii-1,2));
end

% Find gaps below duration threshold
gapIdx = find(gapFrames*stepSize < gapSec*fsHz);

% Set gaps below duration threshold to utterance (logical 1)
if ~isempty(gapIdx)
    for ii = 2:numel(gapIdx)
        out(roi(gapIdx(ii-1),2):roi(gapIdx(ii),1)) = 1;
    end
end


%% APPLY MINIMUM DURATION FOR UTTERANCES
%
%
% Extract utterance lengths in frames
uttFrames = zeros(1,height(roi));
for ii = 1:height(roi)
    uttFrames(ii) = roi(ii,2) - roi(ii,1);
end

% Find utterances below duration threshold
uttIdx = find(uttFrames*stepSize < uttSec*fsHz);

% Set utterances below duration threshold to gaps (logical 0)
if ~isempty(uttIdx)
    for ii = 1:numel(uttIdx)
        out(roi(uttIdx(ii),1):roi(uttIdx(ii),2)) = 0;
    end
end

end