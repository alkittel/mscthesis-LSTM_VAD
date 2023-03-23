function COMB = COMB(x,fsHz,varargin)
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

% First order filter coefficients
beta1 = 10^(-4.8 * P_STFT.stepSize / fsHz);
beta2 = 10^(-2.4 * P_STFT.stepSize / fsHz);
beta3 = 10^(-2.4 * P_STFT.stepSize / fsHz);

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


%% COMPUTE THE MODULATION FEATURES (MPDH, MOD, COMB)
%
%
% Determine number of frames for hold scheme
H = 1 + round(max(0,(P.HSec - P.winSec * P.overlap) / stepSec));

% Dimensionality
nFrames = size(X,2);

% Number of subbands
nBands = size(P.subbands,1);

% Enforce upper limit
H = min(nFrames,H);

% Normalized modulation frequency
w0 = 2 * pi * P.fModHz * P_STFT.stepSize / fsHz;

% Allocate memory
subband = zeros(nBands,nFrames);
MOD = zeros(nBands,nFrames);

% Loop over the number of frequency bands
for kk = 1 : nBands

    % Find frequency bins corresponding to the kk-th band
    idx = fHz >= P.subbands(kk,1) & fHz <= P.subbands(kk,2);

    % Create subband signal
    subband(kk,:) = mean(abs(X(idx,:)));

    % Initialize filter states
    stateIn = 0;
    stateHP = 0;
    stateM  = 0;
    stateNorm = 0;

    % Loop over the number of frames
    for ii = 1 : nFrames

        % High-pass filter
        stateHP = (1 + beta1) * (subband(kk,ii) - stateIn) / 2 + beta1 * stateHP;
        stateIn = subband(kk,ii);

        % Modulation bandpass filter
        stateM = (1 - beta2) * stateHP + beta2 * stateM * exp(1j*w0);

        % Normalization
        if P.bNormMod
            stateNorm = (1 - beta3) * stateHP.^2 + beta3 * stateNorm;

            MOD(kk,ii) = stateM  / sqrt(stateNorm);
        else
            MOD(kk,ii) = stateM;
        end
    end
end

% Modulation-phase difference feature
MPD = -prod(real(MOD),1) + prod(imag(MOD),1);

% Modulation feature
MOD = mean(abs(MOD),1);

% Hold scheme
if P.bCausal

    MPDH = zeros(size(MPD));
    idxH = -H+1:0;
    for ii = 1 : nFrames
        idx = idxH + ii;
        idx = idx(idx > 0);

        MPDH(ii) = max(MPD(idx));
    end
else
    MPDH = ordfilt2(MPD,H,ones(1,H),'symmetric');
end

% Feature combination
COMB = MPDH .* MOD;


end