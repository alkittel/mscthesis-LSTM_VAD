function [bVAD,tSec] = vadGia(x,fsHz,thresholds,P)
%vadGia
% Uses the built-in energy-based detectSpeech function and
% translates output from ROI limits to binary flag array bVAD and corres-
% ponding time vector tSec. Assumes static thresholds have been computed.
%
% INPUTS:
%       x             - input mono signal [1xN]
%       fsHz          - sampling frequency [Hz]
%       thresholds    - see detectSpeech.m
%
%       P - struct containing following detectSpeech parameters:
%             winType       - string defining window type
%             winSize       - frame size in samples
%             overlapLength - window overlap in samples
%             mergeDistance - samples for merging positive speech detection

% OUTPUTS:
%       bVAD - binary array of speech activity flags
%       tSec - time vector
%
% SEE ALSO:
%       detectSpeech.m (for default parameter values)


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Fixed or dynamic thresholds flag
if nargin < 3 || isempty(thresholds)
    dynamicFlag = true;
else
    dynamicFlag = false;
end

% Set default values for parameters
if nargin < 4 || isempty(P)
    P.winType       = 'hann';
    P.winSize       = round(0.03*fsHz);
    P.stepSize      = P.winSize;
    P.mergeDistance = P.winSize*5;
end

% Determine size of input signal
dim = size(x);

% Check if input is mono
if min(dim) > 1
    error('Single-channel input required.')
end


%% DETECT SPEECH
%
%
% Input length
N = numel(x);

% Speech ROI limit indices
if dynamicFlag
    idx = detectSpeech(x,fsHz, ...
        Window=window(P.winType,P.winSize,'periodic'), ...
        OverlapLength=P.winSize-P.stepSize, ...
        MergeDistance=P.mergeDistance);
else
    idx = detectSpeech(x,fsHz, ...
        Thresholds=thresholds, ...
        Window=window(P.winType,P.winSize,'periodic'), ...
        OverlapLength=P.winSize-P.stepSize, ...
        MergeDistance=P.mergeDistance);
end

% Binary array and time vector
bVAD = sigroi2binmask(idx,N);
tSec = (0:N-1).' / fsHz;


end