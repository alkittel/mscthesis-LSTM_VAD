function [bVAD,tSec,P,COMB,gamma] = vadGraf(x,fsHz,varargin)
%vadGraf   Voice activity detector based on modulation-phase differences
%
%USAGE 
%     [bVAD,tSec] = vadGraf(x,fsHz)
%     [bVAD,tSec] = vadGraf(x,fsHz,varargin,...)
%   [bVAD,tSec,P] = vadGraf(x,fsHz,varargin,...)
%
%INPUT ARGUMENTS
%      x : input signal [nSamples x 1]
%   fsHz : sampling frequency in Hertz
% 
%   varargin can be a list of the following parameter/value pairs  
% 
%   =======================================================================
%     PARAMETER | DESCRIPTION                                  | DEFAULT   
%   =======================================================================
%      'winSec' | frame size in seconds                        | 20E-3
%   ------------|----------------------------------------------|-----------
%     'overlap' | window overlap factor                        | 0.5
%   ------------|----------------------------------------------|-----------
%     'initSec' | initial noise-only segment in seconds        | 100E-3
%   ------------|----------------------------------------------|-----------
%        'RSec' | duration across which the LTSV is computed   | 300E-3
%   ------------|----------------------------------------------|-----------
%        'MSec' | duration of Bartlett-Welch averaging         | 100E-3
%   ------------|----------------------------------------------|-----------
%   'bufferSec' | duration across which the LTSV values are    | 1
%               | tracked for minimum and maximum search       |
%   ------------|----------------------------------------------|-----------
%       'alpha' | convex combination parameter                 | 0.3
%   ------------|----------------------------------------------|-----------
%    'fRangeHz' | lower and upper frequency limit of LTSV      | [0 4000]
%   ------------|----------------------------------------------|-----------
%        'hSec' | hangover time constant in seconds            | 0
%   ------------|----------------------------------------------|-----------
%     'bVoting' | binary flag activating VAD voting scheme     | true
%   ------------|----------------------------------------------|-----------
%        'VSec' | duration across which the VAD decisions is   | 300E-3
%               | averaged by the voting scheme                |
%   ------------|----------------------------------------------|-----------
%           'c' | percentage of VAD decisions within 'VSec'    | 0.6
%               | which must be dominated by speech activity   | 
%   ------------|----------------------------------------------|-----------
%      'format' | VAD output format ('samples' or 'frames')    | 'frames'
%   ------------|----------------------------------------------|-----------
%       'bPlot' | binary flag indicating if the VAD output     | false
%               | should be plotted                            |
%   ======================================================================= 
% 
%OUTPUT ARGUMENTS
%   bVAD : binary voice activity decision [nSamples|nFrames x 1]
%   tSec : time axis in seconds [nSamples|nFrames x 1]
%      P : parameter structure
%
%   vadGraf(...) plots the VAD output in a new figure.
% 
%EXAMPLE
%   % Load a chirp signal (y & Fs)
%   load chirp;
% 
%   % Add pink noise at -10 dB and pad the mixture with 1s of noise-only
%   % segments at the beginning and the end
%   mix = genNoisySpeech(y,Fs,'colored_pink',-10,[1 1]);
% 
%   % Perform voice activity detection and plot results 
%   vadGraf(mix,Fs,'initSec',1);
% 
%   See also vadGerven, vadKinnunen, vadMarzinzik, vadRamirez and 
%   vadRamirez2005.
% 
%REFERENCES
%   [1] Ghosh, P. K., Tsiartas, A, and Narayanan, S. (2011). "Robust voice
%       activity detection using long-term signal variability," IEEE
%       Transactions on Audio, Speech, and Language Processing, 19(3),
%       600-613.   

%   Developed with Matlab 9.5.0.944444 (R2018b). Please send bug reports to
%   
%   Author  :  Tobias May, © 2019
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2019/01/14
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
    'initSec'       ...
    'subbands'      ...
    'fModHz'        ...
    'bNormMod'      ...
    'HSec'          ...
    'thres'         ...
    'bufferSec'     ...
    'alpha'         ...
    'bVoting'       ...
    'VSec'          ...
    'c'             ...
    'format'        ...
    'bCausal'       ...
    'bPRC'          ...
    'bPlot'         ...
    'bFlag'         ...
    };

% Parameter settings taken from [1]
paramVal = {'default parameters' mfilename ...
    16E-3       ... % Frame size in seconds
    0.75        ... % Overlap factor
    100E-3      ... % Initial noise-only segment
    [200 2000; 4500 8000] ...
    4           ...
    false       ... % Normalize modulation power
    0.3         ... % Duration of hold scheme
    []          ... % Fixed threshold
    1.25        ... % Duration across which the min & max are tracked
    0.6         ... % Convex combination parameter
    true        ... % Binary flag activating VAD voting scheme
    0.2         ... % Duration across which the VAD decisions are averaged 
    0.5         ... % Percentage of decisions which must be speech-dominated
    'frames'    ... % Output format
    true        ... % Causal flag
    false       ... % PRC threshold
    false       ... % Plot VAD output
    false       ... % Plot VAD output
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

% Determine number of frames for hold scheme
H = 1 + round(max(0,(P.HSec - P.winSec * P.overlap) / stepSec));

% Determine number of frames for buffering maxima and minima
B = 1 + round(max(0,(P.bufferSec - P.winSec * P.overlap) / stepSec));
    
% Calculate initial number of "noise only" frames
nInit = 1 + round(max(0,(P.initSec - P.winSec * P.overlap) / stepSec));


%% CALCULATE THE SPECTRUM
% 
% 
% STFT
[X,framesSec,fHz] = stft(x,fsHz,P_STFT);


%% COMPUTE THE MODULATION FEATURE
% 
% 
% Dimensionality
nFrames = size(X,2);

% Number of subbands
nBands = size(P.subbands,1);

% Enforce upper limit
H = min(nFrames,H);
B = min(nFrames,B);

% Normalized modulation frequency
w0 = 2 * pi * P.fModHz * P_STFT.stepSize / fsHz;

nInit = min(nFrames,nInit);

% Initialize FIFO pointers for buffering
ptrS = 1;
ptrN = 1;

% Allocate memory
thres = 0;
subband = zeros(nBands,nFrames);
MOD = zeros(nBands,nFrames);
gamma = zeros(nFrames,1);
bVAD  = zeros(nFrames,1);
bufferN = zeros(B,1);
bufferS = zeros(B,1);

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


%% POST-PROCESSING
% 
% 
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


%% Adaptive thresholding
% 
%
if ~isempty(P.thres)
    gamma = repmat(P.thres,[1 nFrames]);
    bVAD = COMB > gamma & (1:nFrames) > nInit;
else
    for kk = 1 : nFrames
        
        if kk > nInit && COMB(kk) > thres
            bVAD(kk) = true;
            bufferS(ptrS) = COMB(kk);
            
            % Mover pointer
            ptrS = mod(ptrS,B) + 1;
        else
            bVAD(kk) = false;
            bufferN(ptrN) = COMB(kk);
            
            % Mover pointer
            ptrN = mod(ptrN,B) + 1;
        end
        
        % Compute adaptive threshold (Eq. 7)
        if P.bPRC
            thres = P.alpha * prctile(bufferS,5) + (1 - P.alpha) * prctile(bufferN,95);
        else
            thres = P.alpha * min(bufferS) + (1 - P.alpha) * max(bufferN);
        end
        
        gamma(kk) = thres;
    end
end


%% POST-PROCESS VAD OUTPUT
% 
% 
% Voting scheme
if P.bVoting
    
    % Window size
    V = 1 + round(max(0,(P.VSec - P.winSec * P.overlap) / stepSec));
    V = min(nFrames,V);

    % Indices for voting scheme
    idxV = 0:V-1;
    
    % Loop over the number of frames
    for kk = 1 : nFrames
        
        % Indexing
        idx = idxV + kk;
        idx = idx(idx <= nFrames);
        
        % Voting scheme
        bVAD(kk) = mean(bVAD(idx)) >= P.c;
    end
end


%% RETURN VAD DECISION
% 
% 
% Select output format
switch lower(P.format)
    
    case 'samples'
        
        % Time vector in seconds
        tSec = (1:nSamples).'/fsHz;
        
        % Convert frame-based VAD decision to samples
        bVAD = interp1(framesSec,double(bVAD),tSec,'nearest','extrap');
        
        % Return logical VAD decision
        bVAD = logical(bVAD);
        
    case 'frames'
        
        % Return logical VAD decision
        bVAD = logical(bVAD).';

        % Time vector
        tSec = framesSec.';
        
    otherwise
        error('VAD format "%s" is not supported.',P.format)
end


%% SHOW RESULTS
% 
% 
% If no output is specified
if nargout == 0 || P.bPlot
    
    % Time axis
    t = (1:nSamples)/fsHz;
       
    if isequal(lower(P.format),'frames')
        vIdx = framesSec;
    else
        vIdx = t;
    end
                              
    figure;
    ax(1) = subplot(2,1,1);
    plot(t,x)
    hold on;plot(vIdx,max(abs(x)) * bVAD,'k','LineWidth',2)
    ylabel('Amplitude')
    xlim([0 inf])

    ax(2) = subplot(2,1,2);
    hold on;    
    plot(framesSec,10 * log10(COMB));
    plot(framesSec,10 * log10(gamma),'k--','linewidth',1.5);
    legend({'COMB' '$\gamma$'},'Location','NorthEast',...
        'FontSize',8,'interpreter','latex')
    xlabel('Time (s)')
    ylabel('Amplitude (dB)')
    grid on;
    
    linkaxes(ax,'x');
end


%   ***********************************************************************
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   ***********************************************************************