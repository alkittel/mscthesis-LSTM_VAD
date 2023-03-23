function [x,fsHz] = preprocess(x,fsHz,methods)
%preprocess   Preprocess an audio signal.
%
%USAGE 
%   audio = preprocess(audio,fsHz,method)
%
%INPUT ARGUMENTS
%    audio : audio sigal [nSamples x nChannels]
%     fsHz : sampling frequency in Hertz
%   method : string or cell array defining a list methods
%            '' or 'nothing'  = do nothing
%            'rms'            = RMS normalization
%            'removedc'       = 4th order high-pass with a cutoff frequency
%                               of 20 Hz to remove DC omponents  
%            'lowpass_1000'   = 5th order low-pass with a cutoff frequency 
%                               of 1000 Hz 
%            'lowpass_5000'   = 5th order low-pass with a cutoff frequency 
%                               of 5000 Hz 
%            'whiteningfir'   = first order whitening filter
%            'whiteningfft'   = channel-dependent FFT-based whitening
% 
%OUTPUT ARGUMENTS
%   audio : processed audio signal [nSamples x nChannels]
% 
%NOTE
%   It is possible to combine various preprocessing methods by using a cell
%   array of methods, e.g. method = {'removedc' 'lowpass'}

%   Developed with Matlab 9.1.0.441655 (R2016b). Please send bug reports to
%   
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/01/22
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
if nargin < 3 || isempty(methods); methods = ''; end

% Check if method is a char or a cell array
if ischar(methods) || islogical(methods)
    methods = {methods};
end

% Dimensionality
nMethods = length(methods);


%% PERFORM PRE-PROCESSING
% 
% 
% Loop over number of preprocessing methods
for ii = 1 : nMethods
    
    % Detect underscore to separate method from parameter
    idx = strfind(methods{ii},'_');
    
    if isempty(idx)
        method = methods{ii};
        P = [];
    else
        method = methods{ii}(1:idx-1);
        P = str2double(methods{ii}(idx+1:end));
    end
    
    % Select preprocessing method
    switch lower(method)
        case {'' 'nothing' false}
            % Do nothing ...
            
        case 'downsample'
           % Define new sampling frequency
           if isempty(P)
               P = 2;
           end
           
            % Resample signal
            x = downsample(x,P);

            % Update sampling frequency
            fsHz = round(fsHz / P);
        
        case 'resample'
           % Define new sampling frequency
           if isempty(P)
               P = round(fsHz * 0.25);
           end
           
            % Resample signal
            x = resample(x,P,fsHz);

            % Update sampling frequency
            fsHz = P;
            
        case 'lpc'
           % Define LPC order
           if isempty(P)
               P = round(3/4 * fsHz / 1E3);
           end
           
            % Calculate LPC residual signal
            x = residualLPC(x,fsHz,[],[],[],P);
                
        case 'hp'
           % Define cutoff frequency
           if isempty(P)
               P = 70;
           end
           
            % Create 2nd order high-pass filter
            [b,a] = butter(2, P / (fsHz * 0.5),'high');
            
            % Apply high-pass filter
            x = filter(b,a,x,[],1);
            
        case 'lp'
           % Define cutoff frequency
           if isempty(P)
               P = 5000;
           end
           
            % Create 2nd order low-pass filter
            [b,a] = butter(2, P / (fsHz * 0.5),'low');
            
            % Apply low-pass filter
            x = filter(b,a,x,[],1);

        case 'whiteningfir'
            % Apply 1st order whitening filter
             x = filter([1 -0.97], 1, x);
             
        case 'whiteningfft'
            % Channel-dependent whitening
            for jj = 1 : size(x,2)
                spec        = fft(x(:,jj));
                x(:,jj) = real(ifft(spec./abs(spec)));
            end
            
        otherwise
            error('Pre-processing method "%s" is not supported.',...
                lower(method))
    end
end


%% PLOT AUDIO SIGNAL
% 
% 
% Show signal
if nargout == 0
   figure
   plot(x);
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