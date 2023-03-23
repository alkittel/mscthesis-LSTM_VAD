function [bVAD,tSec,P,LTSV,gamma] = vadGhosh(x,fsHz,varargin)
%vadGhosh   Voice activity detector based on long-term signal variability
%
%USAGE 
%     [bVAD,tSec] = vadGhosh(x,fsHz)
%     [bVAD,tSec] = vadGhosh(x,fsHz,varargin,...)
%   [bVAD,tSec,P] = vadGhosh(x,fsHz,varargin,...)
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
%   vadGhosh(...) plots the VAD output in a new figure.
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
%   vadGhosh(mix,Fs,'initSec',1);
% 
%   See also vadGerven, vadKinnunen, vadMarzinzik, vadRamirez and 
%   vadRamirez2005.
% 
%REFERENCES
%   [1] Ghosh, P. K., Tsiartas, A, and Narayanan, S. (2011). "Robust voice
%       activity detection using long-term signal variability," IEEE
%       Transactions on Audio, Speech, and Language Processing, 19(3),
%       600-613.   

%   Developed with Matlab 9.4.0.813654 (R2018a). Please send bug reports to
%   
%   Author  :  Tobias May, © 2018
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2018/04/25
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
    'bDither'       ...
    'RSec'          ...
    'MSec'          ...
    'bufferSec'     ...
    'alpha'         ...
    'fRangeHz'      ...
    'hSec'          ...
    'bVoting'       ...
    'VSec'          ...
    'c'             ...
    'format'        ...
    'bMT'           ...
    'bPRC'          ...
    'bPlot'         ...
    };

% Parameter settings taken from [1]
paramVal = {'default parameters' mfilename ...
    20E-3       ... % Frame size in seconds
    0.5         ... % Overlap factor
    100E-3      ... % Initial noise-only segment
    true        ... % Binary flag activating dithering
    0.3         ... % Duration across which LTSV should is computed
    0.1         ... % Duration of Bartlett-Welch averaging
    1           ... % Duration across which the LTSV min & max are tracked
    0.3         ... % Convex combination parameter
    [0 4000]    ... % Lower and upper frequency limit of LTSV
    0           ... % Hangover time constant in seconds  
    true        ... % Binary flag activating VAD voting scheme
    0.3         ... % Duration across which the VAD decisions are averaged 
    0.8         ... % Percentage of decisions which must be speech-dominated
    'frames'    ... % Output format
    false       ... % Plot VAD output
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

% Step size in seconds
stepSec = P_STFT.stepSize / fsHz;

% Determine number of frames for LTSV calculation
R = 1 + round(max(0,(P.RSec - P.winSec * P.overlap) / stepSec));

% Determine number of frames for Bartlett-Welch averaging
M = 1 + round(max(0,(P.MSec - P.winSec * P.overlap) / stepSec));

% Determine number of frames for buffering LTSV maxima and minima
B = 1 + round(max(0,(P.bufferSec - P.winSec * P.overlap) / stepSec));
    
% Calculate initial number of "noise only" frames
nInit = 1 + round(max(0,(P.initSec - P.winSec * P.overlap) / stepSec));

% Determine number of frames for hangover scheme
nFramesHO = max(0,1 + round((P.hSec - P.winSec * P.overlap) / stepSec));


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


%% CALCULATE THE POWER SPECTRUM
% 
% 
if P.bMT
    [X,fHz,framesSec] = calcSTFT_MultiTaper(x,fsHz,P.winSec,P.winSec * (1-P.overlap),[],'multitaper_4','analysis');
%     M = 1;
else
    % STFT
    [X,framesSec,fHz] = stft(x,fsHz,P_STFT);
end

% Limit frequency range
X = X(fHz > min(P.fRangeHz) & fHz <= max(P.fRangeHz),:);


%% COMPUTE THE LONG-TERM SIGNAL VARIABILITY (LTSV)
% 
% 
% Dimensionality
nFrames = size(X,2);

% Enforce upper limit
R = min(nFrames,R);
M = min(nFrames,M);
B = min(nFrames,B);

nInit = min(nFrames,nInit);
nFramesHO = min(nFrames,nFramesHO);

% Indices for LTSV calculation
idxLTSV = -R+1:0;

% Indices for Bartlett-Welch averaging
idxAVG = -M+1:0;

% Initialize FIFO pointers for buffering
ptrS = 1;
ptrN = 1;

% Hangover counter
ctrHO = 0;

% Allocate memory
thres = 0;
XPow  = zeros(size(X));
LTSV  = zeros(nFrames,1);
gamma = zeros(nFrames,1);
bVAD  = zeros(nFrames,1);
bufferN = zeros(B,1);
bufferS = zeros(B,1);

% Loop over the number of frames
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
    
    % Apply threshold
    if LTSV(ii) > thres && ii > nInit
        % Speech detected
        bVAD(ii) = true;
        
        % Reset hangover counter
        ctrHO = nFramesHO;
    else
        % Speech pause detected
        bVAD(ii) = false;
    end
        
    % Buffering
    if bVAD(ii)
        % Buffer LTSV metric reflecting speech plus noise
        bufferS(ptrS) = LTSV(ii);
        
        % Mover pointer
        ptrS = mod(ptrS,B) + 1;
    else
        % Buffer LTSV metric reflecting noise 
        bufferN(ptrN) = LTSV(ii);
        
        % Mover pointer
        ptrN = mod(ptrN,B) + 1;
    end
    
    % Compute adaptive threshold (Eq. 7)
    if P.bPRC
        thres = P.alpha * prctile(bufferS,5) + (1 - P.alpha) * prctile(bufferN,95);
    else
        thres = P.alpha * min(bufferS) + (1 - P.alpha) * max(bufferN);
    end
    
    % Hangover scheme
    if bVAD(ii) == 0 && (ctrHO > 0)
        % Prolong speech activity detection
        bVAD(ii) = true;
        
        % Decrease hangover counter
        ctrHO = ctrHO - 1;
    end
    
    % Store adaptive threshold
    gamma(ii) = thres;
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
    for ii = 1 : nFrames
        
        % Indexing
        idx = idxV + ii;
        idx = idx(idx <= nFrames);
        
        % Voting scheme
        bVAD(ii) = mean(bVAD(idx)) >= P.c;
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
    plot(framesSec,10 * log10(LTSV));
    plot(framesSec,10 * log10(gamma),'k--','linewidth',1.5);
    legend({'LTSV' '$\gamma$'},'Location','NorthEast',...
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