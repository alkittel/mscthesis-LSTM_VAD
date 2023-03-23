function signal = zeroPadding(signal,nZeros)
%zeroPadding   Zero-pad input signal along the first dimension
% 
%USAGE
%   signal = zeroPadding(signal,nZeros)
% 
%INPUT ARGUMENTS
%   signal : input signal [nSamples x nChannels]
%   nZeros : number of zeros that should be added to the input signal
%            nZeros = [pre] : pad zeros to the beginning of the signal
%            nZeros = [pre post] : pad zeros to the beginning and the end
%            of the signal 
% 
%OUTPUT ARGUMENTS
%   signal : zero-padded signal [nSampels + sum(nZeros) x nChannels]

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2014
%              Technical University of Denmark (DTU)
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2014/07/17
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 2
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Check if "nZeros" are integers
if any(rem( nZeros,1) ~= 0)
    error('"nZeros" must be integer-valued.')
end

% Check dimensionality 
if numel(size(signal)) > 2
    error('The input signal must be either one or two-dimensional.')
end
    
    
%% PERFORM ZERO-PADDING
% 
% 
% Zero-padding
if any(nZeros > 0)

    % Determine number of channels
    nChannels = size(signal,2);
    
    switch numel(nZeros)
        case 1
            % Preceding zeros
            signal = cat(1,zeros(nZeros,nChannels),signal);
        case 2
            % Preceding and appended zeros
            signal = cat(1,zeros(nZeros(1),nChannels),signal,...
                zeros(nZeros(2),nChannels));
        otherwise
            error('"nZeros" must be either one or two-dimensional.')
    end
end
