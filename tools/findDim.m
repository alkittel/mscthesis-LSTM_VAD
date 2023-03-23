function dim = findDim(x)
%findDim   Find the first non-singleton dimension of a matrix
%   
%USAGE
%   dim = findDim(x)
% 
%INPUT ARGUMENTS
%   x : data matrix
% 
%OUTPUT ARGUMENTS
%   dim : first non-singleton dimension of x

%   Developed with Matlab 8.6.0.267246 (R2015b). Please send bug reports to:
%   
%   Author  :  Tobias May, © 2015
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2015/10/20
%   ***********************************************************************


%% CHECK INPUT ARGUMENTS
% 
% 
% Check for proper input arguments
if nargin < 1 || nargin > 1
    help(mfilename);
    error('Wrong number of input arguments!')
end


%% FIND NON-SINGLETON DIMENSION
% 
% 
% Helper function to find the first non-singleton dimension
dim = find(size(x)~=1,1);
if isempty(dim), dim = 1; end
