function x = istft(X,P,Nx)
%istft   Compute the inverse short-time discrete Fourier transform (ISTFT)
%   A warning message is displayed if the COLA criterion to achieve perfect
%   reconstruction is not satisfied by the chosen STFT parameters.  
% 
%USAGE
%   x = istft(X,P)
%   x = istft(X,P,Nx)
%
%INPUT ARGUMENTS
%    X : single-sided complex STFT matrix consisting of M/2+1 FFT bins, L
%        frames and C channels [M/2+1 x L x C]  
%    P : STFT parameter structure, see configSTFT.m for more details
%   Nx : original signal length in samples of STFT representation X. If
%        specified the output signal x will be trimmed to [Nx x C]
%
%OUTPUT ARGUMENTS
%   x : ouput signal [N + L * R x C] | [Nx x C]
%
%   istft(...) plots the reconstructed output signal in a new figure.
% 
%EXAMPLE
%   % Load signal (y & Fs)
%   load('handel.mat');
% 
%   % Initialize STFT parameter structure 
%   P = configSTFT(Fs);
% 
%   % Compute the STFT
%   Y = stft(y,Fs,P);
% 
%   % Apply some fancy processing in the STFT domain
% 
%   % Reconstruct time-domain signal
%   yHat = istft(Y,P,size(y,1));
%
%   % Samples for which perfect reconstruction is achieved
%   idx = P.offsetCOLA:numel(y)-(P.offsetCOLA-1);
%
%   % Compute the RMS error
%   sqrt(mean(power(y(idx)-yHat(idx),2)))
% 
%   See also stft and configSTFT.

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

% Check if P is a STFT parameter structure
isSTFT(P);

% Check COLA constraint and display a warning message if perfect
% reconstruction is not possible 
if P.bCOLA == false
    warning(['The COLA criterion is not satisfied. To achieve perfect ',...
        'reconstruction, change the sampling frequency, the window ',...
        'length, the window overlap or the window function.'])
end
    

%% INVERSE SHORT-TIME DISCRETE FOURIER TRANSFORM
% 
% 
% Check if nfft size is even
if rem(P.nfft,2) == 0
    % Create double-sided spectrum (do not replicate DC and Nyquist)
    X = cat(1,X,conj(flipud(X(2:end-1,:,:))));
else
    % Create double-sided spectrum (do not replicate DC)
    X = cat(1,X,conj(flipud(X(2:end,:,:))));
end
    
% Get dimensionality
[M,L,C] = size(X);

% Predict total number of samples
T = P.winSize + L * P.stepSize; 

% Select STFT synthesis method
switch(lower(P.ola))
    case 'ola' 
        % =============
        % Reference [1]       
        % =============
        %         
        % Allocate memory
        x = zeros(T,C);

        % Loop over the number of frames
        for ii = 1:L
            % IDFT
            xm = squeeze(real(ifft(X(:,ii,:),P.nfft)));
            
            % Apply scaling factor (a rectangular synthesis window is used)
            xw = xm((1:P.winSize) + P.zeroPadding(1),:) * P.scaling;
            
            % Overlap-and-add individual time frames
            x((1:P.winSize)+(ii-1)*P.stepSize,:) = ...
                x((1:P.winSize)+(ii-1)*P.stepSize,:) + xw;
        end

    case 'wola' 
        % =============
        % Reference [2]       
        % =============
        %                 
        % Allocate memory
        x = zeros(T,C);

        % Loop over the number of frames
        for ii = 1:L
            % IDFT
            xm = squeeze(real(ifft(X(:,ii,:),P.nfft)));
                        
            % Apply synthesis window and scaling factor
            xw = xm((1:P.winSize)+P.zeroPadding(1),:) ...
                .* repmat(P.winS,[1 C]) * P.scaling;
            
            % Overlap-and-add individual time frames
            x((1:P.winSize)+(ii-1)*P.stepSize,:) = ...
                x((1:P.winSize)+(ii-1)*P.stepSize,:) + xw;
        end
        
    case 'allen'
        % =============
        % Reference [1]       
        % =============
        %                 
        % Allocate memory
        x = zeros(T+M,C);
        
        % Loop over the number of frames
        for ii = 1:L
            % IDFT
            xm = squeeze(real(ifft(X(:,ii,:),P.nfft)));
            
            % Apply scaling factor (a rectangular synthesis window is used)
            xw = xm * P.scaling;
            
            % Overlap-and-add individual time frames
            x((1:M)+(ii-1)*P.stepSize,:) = ...
                x((1:M)+(ii-1)*P.stepSize,:) + xw;
        end
        
        % Trim output
        x = x((1:T) + fix((M-P.winSize)/2),:);
    
    case 'mha'
        % =============
        % Reference [3]       
        % =============
        %                 
        % Allocate memory
        x = zeros(T+M,C);
        
        % Loop over the number of frames
        for ii = 1:L
            % IDFT
            xm = squeeze(real(ifft(X(:,ii,:),P.nfft)));
            
            % Apply synthesis window and scaling factor
            xw = xm .* repmat(P.winS,[1 C]) * P.scaling;
            
            % Overlap-and-add individual time frames
            x((1:M)+(ii-1)*P.stepSize,:) = ...
                x((1:M)+(ii-1)*P.stepSize,:) + xw;
        end
        
        % Trim output
        x = x((1:T) + fix((M-P.winSize)/2),:);

    otherwise
        error('Reconstruction method "%s" is not supported.',...
            lower(P.method));
end

% If specified, trim the input x to specified length Nx
if exist('Nx','var')
    if T > Nx
        x(Nx+1:end,:) = [];
    else
        warning(['No trimming is performed because the reconstructed '...
            'signal "x" is shorter than "Nx".'])
    end
end


%% PLOT OUTPUT SIGNAL
%
%
% If no output is specified
if nargout == 0
    n = (0:size(x,1)-1);
    
    figure;
    % Loop over the number of channels
    for ii = 1 : C
        subplot(C,1,ii)
        plot(n,x(:,ii),'color',[0 0.3895 0.9712])
        grid on;
        xlim([n(1) n(end)])
        axis xy
        xlabel('n','interpreter','latex')
        ylabel('Amplitude','interpreter','latex')
    end
end
