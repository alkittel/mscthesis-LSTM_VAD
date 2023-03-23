function STLE = STLE(x,fsHz,varargin)
% ayo. check the graf-ghosh feature function.


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Determine size of input signal
dim = size(x);

% Check if input is mono
if min(dim) > 1
    error('Single-channel input signal required.')
else
    % Number of samples
    nSamples = max(dim);
end


%% CREATE PARAMETER STRUCTURE
%
%
% Parameter names
paramName = {'label' 'fHandle' ...
    'winSec'        ...
    'overlap'       ...
    'winType'       ...
    'subbands'      ...
    'fModHz'        ...
    'bNormMod'      ...
    'HSec'          ...
    'bCausal'       ...
    'RSec'          ...
    'MSec'          ...
    'bDither'       ...
    'fRangeHz'      ...
    };

% Parameter settings taken from [1]
paramVal = {'default parameters' mfilename ...
    15E-3       ... % Frame size in seconds
    0.5         ... % Overlap factor
    'hann'      ... % Window type
    [200 2000; 4500 8000] ...
    4           ...
    false       ... % Normalize modulation power
    0.3         ... % Duration of hold scheme
    true        ... % Causal flag
    0.3         ... % Duration across which LTSV should is computed
    0.1         ... % Duration of Bartlett-Welch averaging
    true        ... % Binary flag activating dithering
    [0 4000]    ... % LTSV frequency limits
    };

% Configure parameter structure
P = configAlgorithm(paramName,paramVal,{'label' 'fHandle'},varargin);


%% INITIALIZE PARAMETERS
%
%
% Configure STFT framework
P_STFT = configSTFT(fsHz,P.winSec,P.overlap,'hann','ola');


%% DITHERING
%
%
% Perform dithering
if P.bDither

    % Assume 16 bit signal representation
    scaling = 2^(16 - 1) - 1;

    % Create dither signal (with triangular PDF)
    dither = rand(nSamples + 1,1) - 0.5;
    dither = dither(1:end-1) + dither(2:end);

    % Dithering and re-scaling
    x = (x * scaling + dither) / scaling;
end


%% CALCULATE THE SHORT TIME LOG ENERGY
%
%
xF = frameData(x,P_STFT.winSize,P_STFT.stepSize,P_STFT.winType);
STLE = log(sum(xF.^2,1));


end