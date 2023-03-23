function [X,t,f,P] = stft(x,fsHz,P)
%stft   Compute the short-time discrete Fourier transform (STFT)
%
%USAGE
%   [X,t,f,P] = stft(x,fsHz)
%   [X,t,f,P] = stft(x,fsHz,P)
%
%INPUT ARGUMENTS
%      x : input signal consisting of N samples and C channels [N x C]  
%   fsHz : sampling frequency in Hertz
%      P : STFT parameter structure, see configSTFT.m for more details
%
%OUTPUT ARGUMENTS
%   X : single-sided complex STFT matrix consisting of M/2+1 FFT bins, L
%       frames and C channels [M/2+1 x L x C] 
%   t : time vector in seconds centered in each frame [1 x L]
%   f : frequency vector in Hertz [M/2+1 x 1]
%   P : STFT parameter structure, see configSTFT.m for more details
% 
%   stft(...) plots the STFT magnitude spectrum in a new figure.
% 
%EXAMPLE 1: Compute the STFT with default parameters
%   % Load signal (y & Fs)
%   load('chirp.mat');
% 
%   % Compute and plot the STFT
%   stft(y,Fs);
% 
%EXAMPLE 2: Compute the STFT with modified parameters
%   % Load signal (y & Fs)
%   load('chirp.mat');
% 
%   % Initialize STFT parameter structure with a window size of 32 ms and 
%   % an overlap factor of 75%
%   P = configSTFT(Fs,32E-3,0.75);
% 
%   % Compute and plot the STFT
%   stft(y,Fs,P);
%
%   See also configSTFT and istft.

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
if nargin < 2 || nargin > 3
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(P); P = configSTFT(fsHz); end

% Check if P is a STFT parameter structure
isSTFT(P);

% Check sampling frequency
if fsHz ~= P.fsHz
    error(['Sampling frequency "fsHz" does not match with the STFT ',...
        'parameter structure.'])
end


%% INITIALIZE PARAMETERS
% 
% 
% Dimensionality
[N,C] = size(x);

% Number of frames
L = max(1,ceil((N-P.overlap)/P.stepSize));   

% Zero-pad input such that it can be divided into an integer number of
% frames 
x = cat(1,x,zeros(round((P.overlap + L * P.stepSize) - N),C));

% Allocate memory for STFT matrix
X = zeros(P.nfft, L, C);     


%% SHORT-TIME DISCRETE FOURIER TRANSFORM
% 
% 
% Loop over the number of frames
for ii = 1:L
    
    % Time segmentation and windowing
    xw = x((1:P.winSize)+(ii-1) * P.stepSize,:) .* repmat(P.winA,[1 C]);
    
    % Zero-padding
    xw = cat(1,zeros(P.zeroPadding(1),C),xw,zeros(P.zeroPadding(2),C));
    
    % Compute and store the DFT
    X(:,ii,:) = fft(xw);
end

% Discard negative frequencies
X = X(1:fix(P.nfft/2)+1,:,:);

% Calculate the time and frequency vectors
t = (P.winSize/2+(0:L-1) * P.stepSize) / fsHz;
f = fsHz/P.nfft * (0:fix(P.nfft/2))';


%% PLOT STFT MAGNITUDE SPECTRUM
%
%
% If no output is specified
if nargout == 0
    
    figure;
    % Loop over the number of channels
    for ii = 1 : C
        subplot(C,1,ii)        
        imagesc(t,f,20*log10(abs(X(:,:,ii))/sum(P.winA.^2)));
        xlim([t(1) t(end)])
        ylim([f(1) fsHz/2])
        colorbar;
        axis xy
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
    end
end
