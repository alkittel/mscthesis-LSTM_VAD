function out = getTestFeatures(x,fsHz,featureMap,varargin)
%getTestFeatures Calculate signal features matching the featureMap for
%testing and evaluating main speaker classifiers
%
%INPUT ARGUMENTS
%              x : input signal [nSamples x 1]
%           fsHz : sampling frequency in Hertz
%     featureMap : map of desired features to extract from x
%
%     varargin can be a list of parameter/value pairs
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
%       Alexander Kittel (student, s164010)
%       Technical University of Denmark
%       Dept. of Health Technology
%
%   Date :
%       2022/12
%   ***********************************************************************

%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Determine size of input signal
dim = size(x);

% Check if input is mono
if min(dim) > 1
    error('Single-channel input signal required.')
end

% This thing apparently
if dim(1,1) == 1
    x = x';
end

desiredFeatures = fields(featureMap);

outCell = cell(numel(desiredFeatures),1);


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
    0.3         ... % Duration across which LTSV should be computed
    0.1         ... % Duration of Bartlett-Welch averaging
    true        ... % Binary flag activating dithering
    [0 4000]    ... % LTSV frequency limits
    };

% Configure parameter structure
P = configAlgorithm(paramName,paramVal,{'label' 'fHandle'},varargin);


%% AUDIO FEATURE EXTRACTOR
%
%
afe = audioFeatureExtractor( ...
    SampleRate = fsHz, ...
    Window = feval(P.winType,P.winSec*fsHz,'periodic'), ...
    OverlapLength = P.overlap*P.winSec*fsHz, ...
    SpectralDescriptorInput = 'melSpectrum');

if isfield(featureMap,'mfcc')
    mfccBands = numel(featureMap.mfcc);
    setExtractorParameters(afe,'melSpectrum',NumBands=mfccBands)
    setExtractorParameters(afe,'mfcc',NumCoeffs=mfccBands)
end

% get desired afe features 
for ii = 1:numel(desiredFeatures)
    if isprop(afe,desiredFeatures{ii})
        afe.(desiredFeatures{ii}) = true;
        outCell{featureMap.(desiredFeatures{ii})(1),:} = extract(afe,x)';
        afe.(desiredFeatures{ii}) = false;
    end
end


%% OTHER FEATURES
%
%
% list of other features considered in this thesis project
otherFeatNames = ["longTermSignalVariability"; ...
    "modulationPhaseDifference"; ...
    "modulationPhaseDifferenceHold"; ...
    "combinedGrafModulationFeature"; ...
    "shortTimeLogEnergy"];

for ii = 1:numel(otherFeatNames)
    if isfield(featureMap,otherFeatNames(ii))
        switch otherFeatNames(ii)
            case 'longTermSignalVariability'
                outCell{featureMap.longTermSignalVariability,:} = ...
                LTSV(x,fsHz,'winSec', P.winSec,'overlap',P.overlap); 
            case 'modulationPhaseDifference'
                outCell{featureMap.modulationPhaseDifference,:} = ...
                    MOD(x,fsHz,'winSec',P.winSec,'overlap',P.overlap);
            case 'modulationPhaseDifferenceHold'
                outCell{featureMap.modulationPhaseDifferenceHold,:} = ...
                    MPDH(x,fsHz,'winSec',P.winSec,'overlap',P.overlap);
            case 'combinedGrafModulationFeature'
                outCell{featureMap.combinedGrafModulationFeature,:} = ...
                    COMB(x,fsHz,'winSec',P.winSec,'overlap',P.overlap);
            case 'shortTimeLogEnergy'
                outCell{featureMap.shortTimeLogEnergy,:} = ...
                    STLE(x,fsHz,'winSec',P.winSec,'overlap',P.overlap);
        end
    end
end


%% GATHER AND OUTPUT
%
%
out = cell2mat(outCell);

end