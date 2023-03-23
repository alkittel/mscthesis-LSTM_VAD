function LTSV = LTSV(x,fsHz,varargin)
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

% Step size in seconds
stepSec = P_STFT.stepSize / fsHz;


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


%% CALCULATE THE SPECTRUM
%
%
% STFT
[X,~,fHz] = stft(x,fsHz,P_STFT);


%% COMPUTE THE LONG-TERM SIGNAL VARIABILITY (LTSV)
%
%
% Dimensionality
nFrames = size(X,2);

% Determine number of frames for LTSV calculation
R = 1 + round(max(0,(0.3 - P.winSec * P.overlap) / stepSec));

% Determine number of frames for Bartlett-Welch averaging
M = 1 + round(max(0,(0.1 - P.winSec * P.overlap) / stepSec));

% Limit frequency range
X = X(fHz > 0 & fHz <= 4000,:);

% Enforce upper limit
R = min(nFrames,R);
M = min(nFrames,M);

% Indices for LTSV calculation
idxLTSV = -R+1:0;

% Indices for Bartlett-Welch averaging
idxAVG = -M+1:0;

% Allocate memory
XPow  = zeros(size(X));
LTSV  = zeros(1,nFrames);

for ii = 1 : nFrames

    % Bartlett-Welch spectral averaging (Eq. 6)
    idxM = idxAVG + ii;
    idxM = idxM(idxM > 0 & idxM <= nFrames);
    XPow(:,ii) = mean(X(:,idxM) .* conj(X(:,idxM)), 2);

    % Calculate R-th order LTSV metric (Eq. 1)
    idxR  = idxLTSV + ii;
    idxR  = idxR(idxR > 0 & idxR <= nFrames);
    sumX  = XPow(:,idxR) ./ repmat(sum(XPow(:,idxR),2),[1 numel(idxR)]);
    zetaK = -sum(sumX .* log(sumX),2);
    LTSV(ii) = var(zetaK);
end

end