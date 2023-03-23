function P = parseCell2Struct(P,inArgs,blocked,bError)
%parseCell2Struct   Parse parameter/value pairs to a structure.
%   The "parameters" are used to find the names of the fields in "P" which
%   should be updated with the corresponding "values".
%
%USAGE
%   P = parseCell2Struct(P,inArgs)
%   P = parseCell2Struct(P,inArgs,blocked,bError)
%      
%INPUT ARGUMENTS
%         P : structure in which field names should be configured according
%             to the parameter/value pairs.  
%    inArgs : cell array consisting of parameter/value pairs (e.g. obtained
%             from MATLAB's varargin) 
%   blocked : list of parameters that should not be updated by the
%             parameter/value pairs (default, blockedFields = '') 
%    bError : if parameter pairs are not recognized, the function reports
%             either an error or a warning message 
%             (default, bError = false) 
% 
%OUTPUT ARGUMENTS
%   P : updated parameter structure
%
%   See also reorganizeStruct.
% 
%EXAMPLES
%   % Create parameter structure
%   P = struct('method','erb','range',[0 4000],'nFilters',32);
% 
%   % Update parameter structure 
%   P = parseCell2Struct(P,{'range',[50 5000]})

%   Developed with Matlab 8.3.0.532 (R2014a). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/02/22 
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 4
    help(mfilename);
    error('Wrong number of input arguments!')
end

% Set default values
if nargin < 2 || isempty(inArgs);  inArgs  = [];    end
if nargin < 3 || isempty(blocked); blocked = '';    end
if nargin < 4 || isempty(bError);  bError  = false; end

% Get number of input arguments to parse
nArg = length(inArgs);

% Check for proper dimension
if rem(nArg,2) ~= 0
    error('Field and parameter arguments must come in pairs.')
end


%% CONFIGURE STRUCTURE
% 
% 
% Loop over number of input arguments
for ii = 1 : 2 : nArg
    % Check if field name should not be changed ...
    if strcmp(blocked,inArgs{ii})
        error(['Field name "',inArgs{ii},'" cannot be changed.'])
    else
        % try
        if isobject(P)
            % Configure P 
            P.(inArgs{ii}) = inArgs{ii+1};
        elseif isfield(P,inArgs{ii})
            P.(inArgs{ii}) = inArgs{ii+1};
        else
            if bError
                % Field does not exist ... return an error
                error('MATLAB:parseCell2Struct',['Field name "',...
                    inArgs{ii},'" is not recognized.'])
            else
                % Field does not exist ... return a warning
                warning('MATLAB:parseCell2Struct',['Field name "',...
                    inArgs{ii},'" is not recognized.'])
            end
        end
    end
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