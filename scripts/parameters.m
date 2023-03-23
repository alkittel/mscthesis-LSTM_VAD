%% SET FIXED PARAMETER VALUES
% 
% 
picpath = '/Users/akittel/Documents/dutten/yr7 thesis/latex/pics/';

% Known signal duration of all recordings
if ~exist('T',"var")
    T = 300;
end

% Sample rate (downsampled)
fsHz = 16E3;

% Signal length
N = fsHz*T; 

% Window type (see valid name list in 'help window')
winType = 'hann';

% Window duration (s) and size (samples)
winSec = 15e-3;
winSize = ceil(winSec * fsHz);

% Window overlap length (samples)
overlapLength = winSize / 2;

% Step size (hop length)
stepSize = winSize - overlapLength;

% Number of frames
nFrames = round((N-overlapLength)/stepSize);

% Minimum utterance and gap durations (s)
utteranceSec = 50E-3;
gapSec = 120E-3;

% Time vectors 
tSec = linspace(0,T,N);
tSecFrame = (winSize/2 + (0:nFrames-1) * stepSize)/fsHz;


%% DISPLAY SET PARAMETERS AND VALUES
% 
% 
disp('WSA parameters:')

disp(['   Recordings are assumed ', num2str(T/60), ' minutes or ', ...
    num2str(N),' samples long.'])

disp(['   Sample rate: ',num2str(fsHz),' Hz'])

disp(['   Window type: ',winType])

disp(['   Window duration: ',num2str(winSec*1E3,2), ...
    ' ms (',num2str(winSize),' samples)'])

disp(['   Window overlap ratio: ', ...
    num2str(overlapLength/winSize,2)])

disp(['   Step size (hop length): ',...
    num2str(stepSize,3)])

disp(['   Number of frames per recorded signal: ',...
    num2str(nFrames,5)])

disp(['   Minimum utterance duration: ',...
    num2str(utteranceSec*1e3,3),' ms']);

disp(['   Minimum inter-utterance pause duration: ',...
    num2str(gapSec*1e3,3),' ms']);