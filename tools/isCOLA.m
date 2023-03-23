function [bCOLA,offset] = isCOLA(winA,winS,R,ola,bPlot)
%isCOLA   Check if STFT windows satisfy the constant overlap-add criterion
%
%USAGE
%   [bCOLA,offset] = isCOLA(winA,winS,R,ola)
%   [bCOLA,offset] = isCOLA(winA,winS,R,ola,bPlot)
%
%INPUT ARGUMENTS
%    winA : STFT analysis window function [nSamples x 1]
%    winS : STFT synthesis window function [nSamples x 1]
%       R : window step size in samples
%     ola : string specifying the overlap-add method for STFT synthesis 
% 
%           'ola' = The overlap-add (OLA) method [1] applies a predefined
%                   analysis window in the analysis stage, while the output
%                   signal is reconstructed by overlapping and adding the
%                   individual windowed output frames using a rectangular
%                   synthesis window.          
% 
%           'wola' = The weighted overlap-add (WOLA) method [2,3] applies a
%                    predefined window twice, once in the analysis and once
%                    in the synthesis stage. The additional synthesis
%                    window aims at reducing discontinuities at the frame
%                    boundaries which can occur due to signal modification
%                    in the STFT domain. Since the window function is
%                    effectively squared, the square-root of the predefined
%                    analysis window is used, such that the multiplication
%                    of the analysis and synthesis window sums up to unity.
%                    Therefore, the WOLA method is restricted to non-
%                    negative window functions.             
% 
%          'allen' = Allen suggested an overlap-add method with symmetric 
%                    zero-padding [1]. The windowed signal is padded 
%                    symmetrically to the next but one higher power
%                    of two, which allows the impulse response of the STFT
%                    manipulation to be as long as the zero-padding without
%                    causing temporal aliasing.         
% 
%            'mha' = The master hearing aid (MHA) employs a weighted
%                    overlap-add method with symmetric zero-padding [4].
%                    This approach is similar to Allen's method, but with a
%                    cosine-tapered synthesis window to smooth out
%                    discontinuities at the frame boundaries which can
%                    occur due to signal manipulation in the STFT domain.    
%
%   bPlot : plot the individual and accumulated analysis/synthesis windows
%           (default, bPlot = false)
%
%OUTPUT ARGUMENTS
%    bCOLA : binary flag indicating if the COLA criterion is satisfied
%   offset : sample offset after which perfect reconstruction is achieved
% 
%   isCOLA(...) plots the individual and accumulated analysis/synthesis
%   windows in a new figure. In addition, a warning message is displayed
%   if the COLA criterion to achieve perfect reconstruction is not
%   satisfied by the chosen STFT parameters.  
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
%   See also isSTFT.

%   Developed with Matlab 9.3.0.713579 (R2017b). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/10/10
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin < 4 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 5 || isempty(bPlot); bPlot = false; end
  

%% PERFORM OLA OR WOLA RECONSTRUCTION
% 
% 
% Determine window length
N = numel(winA);

% Sample offset after which perfect reconstruction is achieved
offset = round((N/R-1) * R + 1);

% Total duration in samples: Transition, 2 windows, Transition
T = round(((N/R-1)*2+2)*R);

% Number of overlapping windows
L = floor((T - (N - R))/R);

% Allocate memory
z = zeros(T,L,3); 

% Loop over the number of frames
for ii = 1 : L
    
    % Window indices
    wIdx = (1:N) + (ii-1) * R;
    
    % Select reconstruction method
    switch(lower(ola))
        case {'ola' 'wola'}
            % ==================
            % References [1,2,3]       
            % ==================
            %                             
            % Normalization constant
            K = R / sum(winA .* winS);

            % Analysis window
            z(wIdx,ii,1) = winA;
            
            % Synthesis window
            z(wIdx,ii,2) = winS;
            
            % Analysis and synthesis window including normalization
            z(wIdx,ii,3) = z(wIdx,ii,3) + (winA .* winS) * K;
            
        case {'allen' 'mha'}
            % ================
            % References [1,4]       
            % ================
            %                             
            % Normalization constant
            K = R / sum(winA);

            % Analysis window
            z(wIdx,ii,1) = winA;
            
            % Synthesis window
            z(wIdx,ii,2) = ones(size(winA));
            
            % Analysis and synthesis window including normalization
            z(wIdx,ii,3) = z(wIdx,ii,3) + (winA .* ones(size(winA))) * K;         
            
        otherwise
            error('Reconstruction method "%s" is not supported.',...
                lower(ola))
    end
end


%% CHECK COLA CRITERION
% 
% 
% Total sum of analysis and synthesis windows
y = sum(z(:,:,3),2);

% RMS error
rmsError = sqrt(mean(power(y(offset:end-offset+1)-1,2)));

% Check if COLA criterion is satisfied
bCOLA = rmsError < 1E-10;


%% PLOT WINDOWS AND COLA CRITERION
%
%
% Plot if specified
if bPlot
    figure;
    hold on;
    h1 = plot(z(:,:,1),'color',[0 0.3895 0.9712],'linestyle','-');
    h2 = plot(z(:,:,2),'color',[0.3725 0.8418 0.6288],'linestyle','--');
    h3 = plot(y,'k-','linewidth',1);
    if bCOLA
        h4 = plot([offset offset; offset T-offset+1; ...
            offset T-offset+1; T-offset+1 T-offset+1]',...
            [0 1;0 0;1 1; 0 1]','--','color',[0 0.5 0],'linewidth',2);
        legend([h1(1) h2(1) h3 h4(1)],{'$w_{A}[n]$' '$w_{S}[n]$' ...
            '$\sum\limits w_{A}[n]w_{S}[n]\kappa$' ...
            'COLA = 1'},'interpreter','latex')
        ylim([0 1])
    else
        legend([h1(1) h2(1) h3],{'$w_{A}[n]$' '$w_{S}[n]$' ...
            '$\sum\limits w_{A}[n]w_{S}[n]$'},'interpreter','latex')
    end
    xlim([1 T])
    xlabel('$n$','interpreter','latex')
    ylabel('Amplitude','interpreter','latex')
    grid on;
end


%% DISPLAY ERROR
% 
% 
% Display error message
if nargout == 0 && ~bCOLA || bPlot && ~bCOLA
    warning(['The COLA criterion is not satisfied. To achieve perfect ',...
        'reconstruction, change the sampling frequency, the window ',...
        'length, the window overlap or the window function.'])
end
