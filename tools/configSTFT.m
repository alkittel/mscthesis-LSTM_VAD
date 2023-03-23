function P = configSTFT(fsHz,winSec,overlap,winType,ola,M,pad,bPow2,bPlot)
%configSTFT   Initialize short-time discrete Fourier transform parameters 
%
%USAGE
%   P = configSTFT(fsHz)
%   P = configSTFT(fsHz,winSec,overlap,winType,ola,M,pad,bPow2,bPlot)
%
%INPUT ARGUMENTS
%      fsHz : sampling frequency in Hertz
%    winSec : window length in seconds, the corresponding window size N in
%             samples is forced to be even: N = 2*round(winSec*fsHz/2) 
%             (default, winSec = 20E-3)  
%   overlap : window overlap factor within the range [0 1), the actual
%             overlap in samples is: O = round(N * overlap) 
%             (default, overlap = 0.75)
%   winType : string specifying the analysis window, of which a periodic
%             version is used: w = window(winType,N,'periodic') 
%                    'rectwin' = rectangular window      
%                       'hann' = hann window        (default)
%                    'hamming' = hamming window 
%                  'blackmann' = blackmann window 
%             'blackmanharris' = 4-term blackman-harris window
%                 'flattopwin' = flat top window
%                 'nuttallwin' = minimum 4-term blackman-harris window
% 
%       ola : string specifying the overlap-add method for STFT synthesis 
%             (default, ola = 'wola')
% 
%               'ola' = The overlap-add (OLA) method [1] applies a
%                       predefined analysis window in the analysis stage,
%                       while the output signal is reconstructed by
%                       overlapping and adding the individual windowed
%                       output frames using a rectangular synthesis window.      
% 
%              'wola' = The weighted overlap-add (WOLA) method [2,3]
%                       applies a predefined window twice, once in the
%                       analysis and once in the synthesis stage. The
%                       additional synthesis window aims at reducing
%                       discontinuities at the frame boundaries which can
%                       occur due to signal modification in the STFT
%                       domain. Since the window function is effectively
%                       squared, the square-root of the predefined analysis
%                       window is used, such that the multiplication of the
%                       analysis and synthesis window sums up to unity.
%                       Therefore, the WOLA method is restricted to
%                       non-negative window functions.           
% 
%             'allen' = Allen suggested an overlap-add method with  
%                       symmetric zero-padding [1]. The windowed signal is
%                       padded symmetrically to the next but one higher
%                       power of two, which allows the impulse response of
%                       the STFT manipulation to be as long as the
%                       zero-padding without causing temporal aliasing.        
% 
%               'mha' = The master hearing aid (MHA) employs a weighted
%                       overlap-add method with symmetric zero-padding [4].
%                       This approach is similar to Allen's method, but
%                       with a cosine-tapered synthesis window to smooth
%                       out discontinuities at the frame boundaries which
%                       can occur due to signal manipulation in the STFT
%                       domain.       
% 
%         M : DFT size, windows are zero-padded if M > N
%             (default, M = 2^ceil(log2(N))   for 'ola' and 'wola'  
%                       M = 2^ceil(log2(N)+1) for 'allen' and 'mha')
%       pad : string specifying type of zero-padding (default, pad = 'sym')
%              'pre' - preceeding zeros
%             'post' - append zeros
%              'sym' - symmetric zero-padding
%     bPow2 : binary flag indicating if the specified window length in
%             seconds should be rounded to the nearest integer power of two
%             (default, bPow2 = false) 
%     bPlot : plot individual and accumulated analysis/synthesis windows
%             (default, bPlot = false)
% 
%OUTPUT ARGUMENTS
%   P : STFT parameter structure used by stft.m and istft.m
%
%   configSTFT(...) plots the individual and accumulated analysis/synthesis
%   windows in a new figure. In addition, a warning message is displayed
%   if the COLA criterion to achieve perfect reconstruction is not
%   satisfied by the chosen STFT parameters. 
% 
%EXAMPLE
%   % Initialize STFT parameter structure with default settings for a
%   % sampling frequency of 16 kHz  
%   configSTFT(16E3)
% 
%REFERENCES
%   [1] Allen, J. B. (1977). "Short term spectral analysis, synthesis and
%       modification by discrete Fourier transform," IEEE Transactions on
%       Acoustic, Speech and Signal Processing, ASSP-25(3), 235-238. 
% 
%   [2] Crochiere, R. (1980). "A weighted overlap-add method of short-time
%       Fourier analysis/synthesis," IEEE Transactions on Acoustic,
%       Speech and Signal Processing, 28(1), 99-102. 
% 
%   [3] Griffin, D. W. and Lim, J. S. (1984). "Signal estimation from
%       modified short-time Fourier transform", IEEE Transactions on 
%       Acoustic, Speech and Signal Processing, ASSP-32(2), 236-243. 
% 
%   [4] Grimm, G., Herzke, T., Berg, D. and Hohmann, V. (2006). "The master
%       hearing aid: A PC-based platform for algorithm development and
%       evaluation", Acta Acustica united with Acustica, 92(4), 618-628.  
% 
%   See also stft and istft.

%   Developed with Matlab 9.2.0.538062 (R2017a). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/09/18
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 1 || nargin > 9
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(winSec);  winSec  = 20E-3;  end
if nargin < 3 || isempty(overlap); overlap = 0.75;   end
if nargin < 4 || isempty(winType); winType = 'hann'; end
if nargin < 5 || isempty(ola);     ola     = 'wola'; end
if nargin < 6 || isempty(M);       M       = [];     end
if nargin < 7 || isempty(pad);     pad     = 'sym';  end
if nargin < 8 || isempty(bPow2);   bPow2   = false;  end
if nargin < 9 || isempty(bPlot);   bPlot   = false;  end

% Supported window functions
winSupported = {'rectwin' 'hann' 'hamming' 'blackman' 'blackmanharris' ...
    'flattopwin' 'nuttallwin'};

% Check selected window
if ~any(strcmp(winSupported,winType))
    error(['The specified window function is not supported. Please ',...
        'select one of the following alternatives: ',...
        repmat('"%s" ',[1 numel(winSupported)])],winSupported{:})
end


%% INITIALIZE PARAMETERS
% 
% 
% Round window size to the nearest integer power of two
if bPow2
    winSec = 2^(round(log2(winSec * fsHz))) / fsHz;
end

% Force the window size to be an even number
N = 2 * round(winSec * fsHz / 2);

% Check overlap factor
if overlap < 0 || overlap >= 1
    error('The overlap factor must be within the range [0 1).')
else
    % Overlap in samples between adjacent windows
    O = round(N * overlap);
end

% Step size between adjacent windows
R = N - O;

% Select STFT method
switch(lower(ola))
    case 'ola' 
        % =============
        % Reference [1]       
        % =============
        %                 
        % Find the next higher power of two
        if ~exist('M','var') || isempty(M)
            M = 2^ceil(log2(N));
        end
        
        % Analysis window function
        if strcmp('rectwin',winType)
            winA = ones(N,1);
        else
            winA = window(winType,N,'periodic');
        end
        
        % Rectangular synthesis window function
        winS = rectwin(N);
        
        % Normalization constant
        K = R / sum(winA .* winS);        
        
    case 'wola' 
        % ===============
        % Reference [2,3]       
        % ===============
        %                 
        % Find the next higher power of two
        if ~exist('M','var') || isempty(M)
            M = 2^ceil(log2(N));
        end
        
        % Analysis window function
        if strcmp('rectwin',winType)
            winA = ones(N,1);
        else
            winA = window(winType,N,'periodic');
        end
        
        % Check if window is non-negative
        if any(winA < 0)
            error(['The WOLA method is restricted to non-negative ',...
                'window functions. Please choose a different window ',...
                'function or change the STFT synthesis method.'])
        end
        
        % Use the square-root window function
        winA = sqrt(winA);
        
        % Square-root synthesis window function
        winS = winA;
        
        % Normalization constant
        K = R / sum(winA .* winS);
        
    case 'allen' 
        % =============
        % Reference [1]       
        % =============
        %                 
        % Find the next but one higher power of two
        if ~exist('M','var') || isempty(M)
            M = 2^ceil(log2(N)+1);
        end
        
        % Analysis window function
        if strcmp('rectwin',winType)
            winA = ones(N,1);
        else
            winA = window(winType,N,'periodic');
        end            
        
        % Rectangular synthesis window function
        winS = rectwin(M);
        
        % Normalization constant
        K = R / sum(winA);
    
    case 'mha' 
        % =============
        % Reference [4]       
        % =============
        %                 
        % Find the next but one higher power of two
        if ~exist('M','var') || isempty(M)
            M = 2^ceil(log2(N)+1);
        end
        
        % Cosine ramps in the synthesis window which correspond to a
        % fraction p of half the zero-padding length (M - N)
        p = 0.25;
       
        % Analysis window function
        if strcmp('rectwin',winType)
            winA = ones(N,1);
        else
            winA = window(winType,N,'periodic');
        end
        
        % Cosine-tapered synthesis window
        Rs   = round((M-N)/2 * p);
        winS = hann(2 * Rs,'periodic');
        winS = [winS(1:Rs); ones(M - 2 * Rs,1); winS(Rs+1:end)];
        
        % Normalization constant
        K = R / sum(winA);
        
    otherwise
        error('STFT synthesis method "%s" is not supported."',lower(ola))
end

% Check if M is larger or equal to N
if M < N
    error('DFT size "M" must be larger or eqal to the window size "N".')
else
    % Zero-padding
    switch(lower(pad))
        case 'pre'
            % Pad preceeding zero
            nZeros = [M - N 0];
        case 'post'
            % Append zeros
            nZeros = [0 M - N];
        case 'sym'
            % Symmetric zero-padding
            nZeros = [floor((M - N)/2) ceil((M - N)/2)];
        otherwise
            error('Zero-padding method "%s" is not supported."',lower(pad))
    end
end


%% CONSTANT OVERLAP ADD (COLA) CRITERION
% 
% 
% Check COLA criterion
if nargout == 0 || bPlot
    % Include plotting
    [bCOLA,offsetCOLA] = isCOLA(winA,winS,R,ola,true);
else
    [bCOLA,offsetCOLA] = isCOLA(winA,winS,R,ola);
end


%% SUMMARIZE STFT PARAMETERS
% 
% 
% Create parameter struct
P = struct('label','STFT parameters','ola',ola,'fsHz',fsHz,'winSize',N,...
    'stepSize',R,'overlap',O,'nfft',M,'zeroPadding',nZeros,...
    'winType',winType,'winA',winA,'winS',winS,'scaling',K,'bCOLA',bCOLA,...
    'offsetCOLA',offsetCOLA);

