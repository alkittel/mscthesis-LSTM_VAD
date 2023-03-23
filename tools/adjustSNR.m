function [mix,target,noise,gain] = adjustSNR(target,noise,snrdB)
%adjustSNR   Adjust the signal-to-noise ratio (SNR) between two signals.
%   The SNR is controlled by adjusting the overall level of the noise. Both
%   signals must be of equal size. 
% 
%USAGE
%   [mix,target,noise,gain] = adjustSNR(target,noise,snrdB)
%
%INPUT ARGUMENTS
%   target : speech signal [nSamples x 1]
%    noise : noise signal  [nSamples x 1]
%    snrdB : SNR in dB
% 
%OUTPUT ARGUMENTS
%      mix : noisy target  [nSamples x 1]
%   target : target signal [nSamples x 1]
%    noise : noise signal  [nSamples x 1]
%     gain : gain factor which is applied to the noise 
% 
%   adjustSNR(...) plots the three signals in a new figure.
% 
%EXAMPLE
%   % Load some music (y & Fs)
%   load handel;
%
%   % Mix music signal with Gaussian noise at 0 dB SNR
%   adjustSNR(y,randn(size(y)),0);
% 
%   See also calcSNR. 

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/10/14
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin ~= 3 
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check dimensions
if min(size(target)) > 1
    error('Single-channel input required.')
end

% Check dimensions
if numel(target) ~= numel(noise)
    error('Speech and noise must be of equal size.')
end


%% ADJUST SNR
% 
% 
% Check if required SNR is finite
if isfinite(snrdB) 
        
    % Compute the current SNR
    snr = calcSNR(target,noise,false,true);
    
    % Compute scaling factor for the noise signal
    gain = sqrt(snr / 10^(snrdB/10));
    
    % Scale the noise to get required SNR
    noise = gain * noise;

elseif isequal(snrdB,inf)
   
    % Handle special case if snrdB = inf, set noise to zero
    gain  = 0;
    noise = noise * gain;
    
elseif isequal(snrdB,-inf)
    
	% Handle special case if snrdB = -inf, set target to zero
    gain   = inf;
    target = target * 0;
    
else
    error(['Specified SNR is not valid: ',num2str(snrdB), ' dB.'])
end
    
% Noisy speech
mix = target + noise;


%% SHOW NOISY TARGET
% 
% 
% If no output is specified
if nargout == 0

    figure;hold on;
    hM = plot(mix);
    hN = plot(noise);
    hS = plot(target);
    grid on;xlim([1 numel(mix)]);
    xlabel('Number of samples')
    ylabel('Amplitude')
    set(hM,'color',[0 0 0])
    set(hS,'color',[0 0.5 0])
    set(hN,'color',[0.5 0.5 0.5])
    legend({'mix' 'noise' 'target'})
end
