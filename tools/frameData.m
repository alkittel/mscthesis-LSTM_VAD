function [frames,win] = frameData(input,winSize,stepSize,win,bPad)
%frameData   Segment input signal into overlapping frames.
% 
%USAGE
%   [frames,win] = frameData(input,winSize);
%   [frames,win] = frameData(input,winSize,stepSize,win,bPad);
% 
%INPUT ARGUMENTS
%      input : input signal [nSamples x 1]
%    winSize : frame size in samples
%   stepSize : step size between adjacent frames in samples 
%              (default, stepSize = round(winSize/2))
%        win : string, function handle or vector defining the window
%              function that is applied to each frame. If win is a string
%              or a function handle, a periodic window function will be
%              used (default, win = 'rectwin') 
%       bPad : perform zero-padding to make the input signal evenly
%              divisible by an integer number of frames 
%              (default, bPad = true)
% 
%OUTPUT ARGUMENTS
%   frames : 2D matrix of overlapping frames where each column represents
%            one frame [winSize x nFrames]  
%      win : window function [winSize x 1]
% 
%EXAMPLE
%   % Input signal 
%   input = randn(1E3,1); 
% 
%   % Segment input into frames of 256 samples with 128 samples overlap
%   frames = frameData(input,256);

%   Developed with Matlab 7.9.0.529 (R2009b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2009-2015 
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2009/11/23
%   v.0.2   2010/05/18 automatically determine the number of frames
%   v.0.3   2013/08/19 added zero-padding
%   v.0.4   2015/02/20 added documentation
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
if nargin < 3 || isempty(stepSize); stepSize = round(winSize / 2); end
if nargin < 4 || isempty(win);      win      = 'rectwin';          end
if nargin < 5 || isempty(bPad);     bPad     = true;               end
    
% Determine size of input signal
dim = size(input);
    
% Check if input is mono
if min(dim) > 1
    error('Single-channel input signal is required.')
else
    % Ensure colum vector
    input = input(:);
    
    % Number of samples
    nSamples = max(dim);
end

% Check if window is a vector in samples or a string | function handle
if ischar(win) || isa(win,'function_handle')
    % Check if a window function other than the rectangular should be used
    if isequal(win,'rectwin') || isequal(win,@rectwin)
        win = ones(winSize,1);
        bWindowing = false;
    else
        win = window(win,winSize,'periodic');
        bWindowing = true;
    end
else
    if winSize ~= length(win)    
       error('Mismatch between "winSize" and length of "win".') 
    else
        bWindowing = true;
    end
end


%% ZERO-PADDING
% 
% 
% Overlap
overlap = winSize-stepSize;

% Compute the number of frames
nFrames = (nSamples-overlap)/(stepSize);

% Append zeros
if bPad
    % Number of frames (at least one)
    nFrames = max(1,ceil(nFrames));
    
    % Compute the number of required zeros
    nZeros = (stepSize * nFrames + overlap) - nSamples;
    
    % Pad input signal with zeros
    input = [input; zeros(nZeros,1)];
else
    % No zero padding, exclude elements that do not fit into last frame
    nFrames = max(0,floor(nFrames));
end


%% FRAME INPUT SIGNAL
% 
% 
% Frame indices
frameIdx = repmat((1:winSize)',[1 nFrames]) + ...
    repmat((0:nFrames-1) * stepSize,[winSize 1]);

% Perform framing
frames = input(frameIdx);


%% WINDOWING
% 
% 
% Apply window function
if bWindowing
    frames = bsxfun(@times,frames,win);
end
