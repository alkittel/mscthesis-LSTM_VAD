function [feats,fullFeatureMap] = graf_ghosh_features(x,fsHz,varargin)
%graf_ghosh_features Extract features from VAD techniques described in [1]
% and [2].
%
%USAGE
%     [feats,featureMap] = graf_ghosh_features(x,fsHz)
%     [feats,featureMap] = graf_ghosh_features(x,fsHz,varargin,...)
%
%INPUT ARGUMENTS
%      x : input signal [nSamples x 1]
%   fsHz : sampling frequency in Hertz
%
%   varargin can be a list of parameter/value pairs
%
%REFERENCES
%   [1] Ghosh, P. K., Tsiartas, A, and Narayanan, S. (2011). "Robust voice
%       activity detection using long-term signal variability," IEEE
%       Transactions on Audio, Speech, and Language Processing, 19(3),
%       600-613.
%
%   [2] Graf, S., Herbig, T., Buck, M., & Schmidt, G. (2016). "Voice
%       activity detection based on modulation-phase differences."
%       Speech Communication - 12. ITG-Fachtagung Sprachkommunikation,
%       80–84.
%
%   Developed with Matlab 9.12.0.2039608 (R2022a)
%
%   Authors :
%       Tobias May, © 2018
%       Technical University of Denmark
%       tobmay@elektro.dtu.dk
%
%       Alexander Kittel (student, s164010)
%       Technical University of Denmark
%       Dept. of Health Technology
%
%   Date :
%       2022/11/17
%   ***********************************************************************


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


%% CALCULATE THE SHORT TIME LOG ENERGY
%
%
xF = frameData(x,P_STFT.winSize,P_STFT.stepSize,P_STFT.winType);
STLE = log(sum(xF.^2,1));


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


%% COMPUTE THE LONG-TERM SIGNAL VARIABILITY (LTSV)
%
%
%
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


%% OUTPUT FEATURES
%
%
% Combine features
feats = [LTSV;MOD;MPDH;COMB;STLE];
fullFeatureMap = struct( ...
    'longTermSignalVariability',1, ...
    'modulationPhaseDifference',2, ...
    'modulationPhaseDifferenceHold',3, ...
    'combinedGrafModulationFeature',4, ...
    'shortTimeLogEnergy',5);


end