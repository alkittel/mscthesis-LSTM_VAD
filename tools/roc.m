function R = roc(target,test,dim)
%roc   Calculate receiver operating characteristic (ROC) statistics
%
%USAGE 
%   R = roc(target,test)
%   R = roc(target,test,dim)
%
%INPUT ARGUMENTS
%   target : binary target pattern 
%     test : binary test pattern 
%      dim : dimension along which the ROC statistics should be computed
%            (default, dim = findDim(target))
% 
%OUTPUT ARGUMENTS
%   R : ROC structure
%       .TPR  - true positive (hit) rate               
%       .FNR  - false negative (miss) rate              
%       .FPR  - false positive (false alarm) rate              
%       .TNR  - true negative (correct rejection) rate               
%       .P    - precision                        
%       .ACC  - accuracy    
%       .BACC - balanced accuracy    
%       .F    - F-measure                        
%       .HFA  - hit rate - false alarm rate      
%       .FAR  - overall false alarm error norm
%       .MCC  - Matthews correlation coefficient 
% 
%REFERENCES
%   [1] Fawcett, T. (2006). "An introduction to ROC analysis," Pattern
%       Recognition Letters, 27(8), 867-874.   

%   Developed with Matlab 9.4.0.813654 (R2018a). Please send bug reports to
%   
%   Author  :  Tobias May, © 2018
%              Technical University of Denmark
%              tobmay@elektro.dtu.dk
%
%   History :
%   v.0.1   2018/04/30
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
if nargin < 3 || isempty(dim); dim = findDim(target); end

% Check dimensionality
if size(target) ~= size(test)
    error('Target and test pattern must be of equal size.')
end

% Check if target and test are binary
if mean((target(:) == 0) + (target(:) == 1)) < 1 || ...
        mean((test(:) == 0) + (test(:) == 1)) < 1
    error('Target and test pattern must be binary.')
end


%% ROC ANALYSIS
% 
% 
% Detect activity in target and test matrix
bTarget = target == 1;
bTest   = test == 1;

% ROC statistics
TP = sum( bTest &  bTarget, dim); % True positive
FP = sum( bTest & ~bTarget, dim); % False positive
TN = sum(~bTest & ~bTarget, dim); % True negative
FN = sum(~bTest &  bTarget, dim); % False negative

POS = TP  + FN;  % Total positives (w.r.t. target)
NEG = FP  + TN;  % Total negatives (w.r.t. target)
N   = POS + NEG; % Total items

% Calculate rates
R.TPR = TP ./ POS;
R.FNR = FN ./ POS;
R.FPR = FP ./ NEG;
R.TNR = TN ./ NEG;


%% CALCULATE METRICS
% 
% 
% Precision
R.P = TP ./ (TP + FP);

% Accuracy
R.ACC = (TP + TN) ./ N;

% Balanced accuracy
R.BACC = 0.5 * (R.TPR + R.TNR);

% The F-measure is the harmonic mean of precision and recall
R.F = 2 .* TP ./ (2 .* TP + FP + FN);

% Hit rate - false alarm rate
R.HFA = R.TPR - R.FPR;

% Overall false alarm error norm:
R.FAR = sqrt((1 - R.TPR).^2 + (1 - R.TNR).^2);

% Matthews correlation coefficient
R.MCC = (TP .* TN - FP .* FN) ./ ...
    sqrt( (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN) );
