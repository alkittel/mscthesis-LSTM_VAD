function P = configAlgorithm(name,val,blocked,parse)
%configAlgorithm   Initialize default parameters and customize settings
% 
%USAGE
%   P = configAlgorithm(name,val)
%   P = configAlgorithm(name,val,blocked,parse)
% 
%INPUT ARGUMENTS
%      name : cell array of parameter names
%       val : cell array of default parameter values
%   blocked : list of parameter names that should not be updated by the
%             parameter/value pairs (default, blocked = '')
%     parse : can be one of the following options:
%             1) cell array of parameter/value pairs that should be
%                customized (e.g. obtained from MATLAB's varargin) 
%             2) parameter structure from a previous function call
%                containing all parameter names.
% 
%OUTPUT ARGUMENTS
%   P : parameter structure
% 
%   See also parseCell2Struct.  
% 
%EXAMPLES
%   % Create cell array of parameter names
%   name = {'method' 'range' 'nFilters'};
% 
%   % Create cell array of default parameter values
%   val = {'erb' [0 4000] 32};
% 
%   % Create default parameter structure
%   P1 = configAlgorithm(name,val)
% 
%   % Create default parameter structure with customized settings
%   P2 = configAlgorithm(name,val,[],{'range',[50 5000]})

%   Developed with Matlab 9.0.0.341360 (R2016a). Please send bug reports to
%   
%   Author  :  Tobias May, © 2016
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2016/05/19
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 2 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 3 || isempty(blocked); blocked = ''; end
if nargin < 4 || isempty(parse);   parse   = {}; end

% Check for proper inputs
if ~iscell(name)
    error('"name" must be a cell array.')
end
if ~iscell(val)
    error('"val" must be a cell array.')
end
if numel(name) ~= numel(val)
    error('"name" and "val" must be of equal size.')
end


%% CHECK IF PARSE CONTAINS A STRUCT
% 
% 
% Verify if structure contains all parameter names 
if numel(parse) == 1 && isstruct(parse{1})
    
    % Allocate memory
    bExist = false(numel(name),1);
    
    % Loop over all parameter names
    for ii = 1 : numel(name)
        if isfield(parse{1},name{ii})
            bExist(ii) = true;
        else
            if any(strcmp(blocked,name{ii}))
                bExist(ii) = true;
            end
        end
    end
    
    % Report error if not all parameter names were detected
    if any(~bExist)
        error(['The following parameters are missing in the input ',...
            'structure:',repmat(' "%s"',[1 sum(~bExist)])],...
            name{~bExist})
    else
        P = parse{1};
        return;
    end
end


%% CREATE PARAMETER STRUCT
% 
% 
% Create parameter struct with default settings
try
    P = cell2struct(val,name,2);
catch ME
    error(['Parameter configuration failed: %s Please review ',...
        'the list of controllable parameters.'],ME.message);
end

% Parse additional input arguments (return error of field does not exist!)
try
    P = parseCell2Struct(P,parse,blocked,true);
catch ME
    if isfield(P,'fHandle')
        error(['Parameter configuration of "%s" failed: %s Please review ',...
            'the list of controllable parameters by typing the following ',...
            'in the command window: help %s '],P.fHandle,ME.message,P.fHandle);
    else
        error(['Parameter configuration failed: %s Please review ',...
            'the list of controllable parameters.'],ME.message);
    end
end

% Add comment if the field label exist
if ~isempty(parse) && isfield(P,'label')
    P.label = [P.label ' + varargin'];
end