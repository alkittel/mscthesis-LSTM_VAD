function [bSTFT,fMissing] = isSTFT(P)
%isSTFT   Check if input is a STFT parameter structure
%
%USAGE
%   bSTFT = isSTFT(P)
%
%INPUT ARGUMENTS
%   P : STFT parameter structure 
% 
%OUTPUT ARGUMENTS
%      bSTFT : binary flag indicating if P is a STFT parameter structure
%   fMissing : cell array with missing field names
% 
%   isSTFT(...) will display an error message and the missing field names
%   if P is not a STFT parameter structure.    
% 
%   See also isCOLA.

%   Developed with Matlab 9.3.0.713579 (R2017b). Please send bug reports to
%
%   Author  :  Tobias May, © 2017
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2017/10/24
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
%
%
% Check for proper input arguments
if nargin ~= 1
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% COMPARE STRUCT FIELDS
% 
% 
% Required field names
fieldNames = {'label','ola','fsHz','winSize','stepSize','overlap',...
    'nfft','zeroPadding','winType','winA','winS','scaling','bCOLA',...
    'offsetCOLA'};

% Initialize flags 
bField = false(numel(fieldNames),1);
bSTFT = false;

% Create if P is a struct
if isstruct(P)
	% Get all field names
    fields = fieldnames(P);
    
    % Loop over the number of fields
    for ii = 1 : numel(fields)
        bField(strcmpi(fieldNames,fields{ii})) = true;
    end
    
    % Set STFT flag if all fields exist
    if all(bField)
       bSTFT = true;
    end
end

% Return missing field names
fMissing = fieldNames(~bField);

% Display error message
if nargout == 0 && ~bSTFT
    error(['"P" is not a STFT parameter structure. The following ',...
        'fields are missing: ',repmat('''%s'' ',[1 numel(fMissing)])],...
        fMissing{:})
end
