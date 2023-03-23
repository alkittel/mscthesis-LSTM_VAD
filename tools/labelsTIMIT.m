function [marker,labels,vad,labelsVAD] = labelsTIMIT(fileName)
%labelsTIMIT   Create label vector for TIMIT sentences

% Initialize persistent memory
persistent PERphonetics PERlabels

% Initialize phonetic alphabet
if isempty(PERphonetics) || isempty(PERlabels)

    % Allocate cell arrays
    phonetics = cell(4,1);
    symbols   = cell(8,1);
    labels    = cell(1,1);
    voicing   = cell(1,1);
    
    % Others
%     symbols{1} = {'pau' 'epi' 'h#' '1' '2'};
    symbols{1} = {'pau' 'epi' 'h#'};
    labels{1}  = 'others';
    voicing{1} = 'silent';
    
    % Stops
    symbols{2} = {'b' 'd' 'g' 'p' 't' 'k' 'dx' 'q'};
    labels{2}  = 'stops';
    voicing{2} = 'unvoiced';
    
    % Closure of stops
    symbols{3} = {'bcl' 'dcl' 'gcl' 'pcl' 'tck' 'kcl' 'tcl'};
    labels{3}  = 'closure';
    voicing{3} = 'unvoiced';
    
    % Affricatives
    symbols{4} = {'jh' 'ch'};
    labels{4}  = 'affricates';
    voicing{4} = 'unvoiced';
    
    % Fricatives
    symbols{5} = {'s' 'sh' 'z' 'zh' 'f' 'th' 'v' 'dh'};
    labels{5}  = 'fricatives';
    voicing{5} = 'unvoiced';
    
    % Nasals
    symbols{6} = {'m' 'n' 'ng' 'em' 'en' 'eng' 'nx'};
    labels{6}  = 'nasals';
    voicing{6} = 'voiced';
    
    % Semivowels
    symbols{7} = {'l' 'r' 'w' 'y' 'hh' 'hv' 'el'};
    labels{7}  = 'semivowels';
    voicing{7} = 'voiced';
    
    % Vowels
    symbols{8} = {'iy' 'ih' 'eh' 'ey' 'ae' 'aa' 'aw' 'ay' 'ah' 'ao' 'oy' 'ow' 'uh' 'uw' 'ux' 'er' 'ax' 'ix' 'axr' 'ax-h'};
    labels{8}  = 'vowels';
    voicing{8} = 'voiced';
    
    % Create phonetic alphabet
    for ii = 1 : numel(symbols)
        phonetics{1} = cat(2,phonetics{1},symbols{ii});
        phonetics{2} = cat(2,phonetics{2},repmat(labels(ii),[1 numel(symbols{ii})]));
        phonetics{3} = cat(2,phonetics{3},ii * ones(1,numel(symbols{ii})));
        phonetics{4} = cat(2,phonetics{4},repmat(voicing(ii),[1 numel(symbols{ii})]));
    end
    
    % Update persistent variable
    PERphonetics = phonetics;
    PERlabels    = labels;
    
    % Clean up
    clear symbols voicing ii;
else
    % Use persistent variable
    phonetics = PERphonetics;
    labels    = PERlabels;
end


%% READ PHONETICS
% 
% 
% Get TIMIT labels and temporal marker
phn = readPhonetics(fileName);

% Number of samples in the phn file
nSamples = phn{end,1}(end);

% Number of phonemes
nPhonemes = size(phn,1);

% Allocate memory
marker = zeros(nSamples,1);
vad    = zeros(nSamples,1);

% VAD labels
labelsVAD = [-1 0 1];
strVAD    = {'silent' 'unvoiced' 'voiced'};


%% CREATE SAMPLE-BASED TIMIT LABELS 
% 
% 
% Loop over number of phonemes
for ii = 1:nPhonemes
    
    % Detect phoneme
    bPhon = strcmpi(deblank(phn{ii,2}),phonetics{1});
    
    if any(bPhon)
        
        % Current class
        currClass = phonetics{3}(bPhon);
        
        % Current time segment
        currSeg = phn{ii,1}(1):phn{ii,1}(2);
        
        % Construct labels
        marker(currSeg) = currClass * ones(numel(currSeg),1);
        
        % Current VAD class
        currVADIdx = strcmp(phonetics{4}{bPhon},strVAD);

        % VAD
        if any(currVADIdx)
            vad(currSeg) = labelsVAD(currVADIdx) * ones(numel(currSeg),1);
        else
            error('VAD class ''%s'' is not recognized!',phonetics{4}(bPhon))
        end
    else
        error('Phoneme ''%s'' is not recognized.',deblank(phn(ii,:)))
    end
end


%% CHECK IF PHONETIC INFORMATION MATCHES THE DURATION OF THE WAVE FILE
% 
% 
% Get details on wav file
dim = readAudio(fileName,'info');

% Difference between wave file and phn signal
nSamplesDiff = dim(1) - nSamples;

% Pad marker, in case the last information is "silent"
if nSamplesDiff > 0

    if strcmp('silent',PERphonetics{4}(find(PERphonetics{3}==marker(end),1)))
        marker = cat(1,marker,repmat(marker(end),[nSamplesDiff 1]));
        vad = cat(1,vad,repmat(vad(end),[nSamplesDiff 1]));
    else
       error(['Cannot replicate a marker different from "silent" to',...
           ' match duration between wave file and phonetic labels.'])
    end
end


function phn = readPhonetics(fileName)
%readPhonetics   Read phonetic labels from TIMIT database.

% Get filename
[root,name,suffix] = fileparts(fileName); %#ok

% Add phn suffix
namePhn = [root,filesep,name,'.phn'];

if exist(namePhn,'file')
    fidw = fopen(namePhn,'r');
    
    phn = cell(0,0);
    
    if fidw > 0
        while 1
            tline = fgetl(fidw); % read an input line
            if ~ischar(tline)
                break
            end
            [wtim, ntim, ee, nix] = sscanf(tline,'%d%d',2); %#ok
            if ntim==2
                phn{end+1,1} = max(wtim(:)',1); %#ok, prevent 0
                phn{end,2}   = strtrim(tline(nix:end));
            end
        end
        fclose(fidw);
    end
else
   error('Phonetical file ''%s'' does not exist.',namePhn) 
end