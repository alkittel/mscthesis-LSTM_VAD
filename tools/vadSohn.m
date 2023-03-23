function [pVAD,tSec] = vadSohn(x,fsHz,VAD,P)
%vadSohn
%
% INPUTS:
%       x          - input mono signal [1xN]
%       fsHz       - sampling frequency [Hz]
%       VAD        - voiceActivityDetector system object
%       P          - struct containing following parameters:
%                     winType  - string defining window type
%                     winSize  - frame size in samples
%                     stepSize - window overlap in samples%
%
% OUTPUTS:
%       pVAD - probability of speech activity per frame
%       tSec - time vector 
%
% SEE ALSO:
%       voiceActivityDetector documentation


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values for parameters
if nargin < 5 || isempty(P)
    P.winType  = 'hann';
    P.winSize  = round(0.03*fsHz);
    P.stepSize = P.winSize; % no overlap
end

% Determine size of input signal
dim = size(x);

% Check if input is mono
if min(dim) > 1
    error('Single-channel input required.')
end

% Frame input signal
framedInput = frameData(x,P.winSize,P.stepSize, ...
    P.winType,true);


%% DETECT SPEECH
% 
% 
% Speech presence probability of each signal frame
pVAD = zeros(size(framedInput,2),1); % prealloc
for ii = 1:numel(pVAD)
    pVAD(ii) = VAD(framedInput(:,ii));
end

% Time vector
nFrames = size(framedInput,2);
tSec = (P.winSize/2 + (0:nFrames-1).' * P.stepSize)/fsHz;

end