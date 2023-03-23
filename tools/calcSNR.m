function snr = calcSNR(target,noise,bdB,bRemoveDC,dim)
%calcSNR   Calculate the signal-to-noise ratio (SNR) between two signals.
%   The DC component of both signals can be removed prior to calculating
%   the SNR by using the flag 'bRemoveDC'.   
% 
%USAGE
%   snr = calcSNR(target,noise)
%   snr = calcSNR(target,noise,bdB,bRemoveDC,dim)
%
%INPUT ARGUMENTS
%      target : target signal [nSamples x nChannels]
%       noise : noise signal [nSamples x nChannels]
%         bdB : binary flag indicating if the SNR should be specified in
%               dB. If false, the linear SNR is computed
%               (default, bdB = true)
%   bRemoveDC : binary flag indicating if the DC component of both signals
%               should be removed prior to SNR calculation 
%               (default, bRemoveDC = true)
%         dim : dimension across which the SNR should be computed. By
%               default, the SNR is calculated across the first
%               non-singleton dimension of the target signal 
%               (default, dim = findDim(target))   
% 
%OUTPUT ARGUMENTS
%   snr : SNR (either linear or in dB) [nSamples x 1] | [1 x nChannels]
% 
%   See also adjustSNR.

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/20
%   v.0.2   2015/02/21 added multi-channel support
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS  
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 5
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default parameter
if nargin < 3 || isempty(bdB);       bdB       = true;            end
if nargin < 4 || isempty(bRemoveDC); bRemoveDC = true;            end
if nargin < 5 || isempty(dim);       dim       = findDim(target); end

% Check dimensions
if size(target) ~= size(noise)
    error('"Target" and "noise" must be of equal size.')
end


%% CALCULATE THE SNR
% 
% 
% DC removal
if bRemoveDC
    target = bsxfun(@minus,target,mean(target,dim));
    noise  = bsxfun(@minus,noise, mean(noise ,dim));
end
    
% Energy of target and noise signal 
energyTarget = sum( target.^2 , dim );
energyNoise  = sum( noise .^2 , dim );

% Compute the SNR 
if bdB
    % dB domain
    snr = 10 * log10(energyTarget ./ energyNoise);
else
    % Linear domain
    snr = energyTarget ./ energyNoise;
end
