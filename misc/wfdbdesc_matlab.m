function [ rec_info, sig_info, comment ] = wfdbdesc_matlab(fn,datapath)
% Implements WFDBDESC using only MATLAB textscan functions - no java.

% rec_info 
%       recname, numseg, numsig,...
%       fs, fc, fc_base,...
%       numsamp, basetime, basedate

% sig_info has columns: 
%       gain, baseline, unit,...
%       adcresolution, adczero, initialvalue,...
%       checksum, blocksize, desc

%% load in the .hea file
fp = fopen([datapath fn '.hea']);
info = textscan(fp,'%s','delimiter','\n');
info = info{1};
fclose(fp);

%% move comments to their own variable
idxComment = strncmp(info,'#',1);
comment = info(idxComment);
info = info(~idxComment);

%% get info about the record

%=== default values
numseg=1;
numsig=0;
fs = 250; % default frequency, DEFREQ, defined in <wfdb/wfdb.h>
fc = fs;
fc_base = 0;
numsamp = 0;
basetime = '00:00:00';
basedate = [];

recsplit = regexp(info{1},' ','split');
S = numel(recsplit);
recname = recsplit{1};
idxSlash = regexp(recname,'/','once');
if ~isempty(idxSlash) % extract the number of segments in the record
    numseg = str2double(recname(idxSlash(1)+1:end));
    recname = recname(1:idxSlash(1)-1);
end

if S>=2
    numsig = str2double(recsplit{2});
end
if S>=3
    % check for optional counter frequency delimited by a slash
    idxSlash = regexp(recsplit{3},'/','once');
    if isempty(idxSlash) % no counter frequency, only sampling freq
        fs = sscanf(recsplit{3},'%f');
        fc = fs;
        fc_base = 0;
    else % also parse counter frequency
        fs = sscanf(recsplit{3}(1:idxSlash(1)-1),'%f'); % sampling frequency
        fc = recsplit{3}(idxSlash(1)+1:end); % counter frequency
        
        % check for optional base counter value
        fc_base = regexp(fc,'\([0-9]*\)','match');
        if ~isempty(fc_base)
            fc_base = str2double(fc_base{1}(2:end-1)); % removes the parentheses
            fc = sscanf(fc(1:regexp(fc,'\(','once')-1),'%f');
        else
            fc = sscanf(fc,'%f');
        end
    end
end
if S>=4
    numsamp = str2double(recsplit{4});
end
if S>=5
    basetime = recsplit{5};
end
if S>=6
    basedate = recsplit{6};
end

rec_info = {recname, numseg, numsig,...
    fs, fc, fc_base,...
    numsamp, basetime, basedate};

%% get info about each individual signal
header = cell(1,numel(info)-1);
sig_info = cell(numel(info)-1,9);
for k=2:numel(info)
    tmp = info{k}(end:-1:1); % flip the string to get the last segment
    idxSpace = regexp(tmp,' ','once');
    header{k} = tmp(1:idxSpace(1));
    % flip the header back to the correct order
    header{k} = header{k}(end:-1:1);
    
    % get the original string again
    tmp = info{k};
    
    %% default parameters according to HEADER(5) definition
    
    % default the rest of the parameters
    gain = 200; % default gain defined in wfdb/wfdb.h
    baseline = 0;
    unit = [];
    adcresolution=12;
    adczero = 0;
    initialvalue = 0;
    checksum = 0;
    blocksize = 0;
    desc = [];
    
    %% now start scanning through a signal line to interpret the input
    % this section loops through whitespace delimited fields
    % as soon as whitespaces stop, it defaults the remaining parameters
    % there must be *at least* two columns, the rest are optional, but the
    % order of the optional arguments is fixed
    
    sigsplit = regexp(tmp,' ','split');
    S = numel(sigsplit);
    
    % first get the signal name
    signame = sigsplit{1};
    
    % next, get the format
    frmt = sigsplit{2};    
    
    if S>=3
        % if s>2, we have a gain
        gain = sigsplit{3};
        
        % the above can contain baseline/ADC resolution
        
        baseline = regexp(gain,'\([0-9]*\)','match');
        if ~isempty(baseline)
            baseline = str2double(baseline{1}(2:end-1)); % removes the parentheses
        else
            baseline = 0;
        end
        idxSlash = regexp(gain,'/','once');
        if ~isempty(idxSlash)
            unit = gain(idxSlash(1)+1:end);
            gain = str2double(gain(1:idxSlash(1)-1));
        else
            unit = [];
            gain = str2double(gain);
        end
    end
    if S>=4
        adcresolution = str2double(sigsplit{4});
        if adcresolution==0; adcresolution=12; end
    end
    if S>=5
        adczero = str2double(sigsplit{5});
        if S==5
            initialvalue = adczero; % special case 4 fields, initial value is equal to ADC zero value
        end
    end
    if S>=6
        initialvalue = str2double(sigsplit{6});
    end
    if S>=7
        checksum = str2double(sigsplit{7});
    end
    if S>=8
        blocksize = str2double(sigsplit{8});
    end
    %== 8 fields ==%
    if S>=9
        desc = sigsplit{9};
    end
    
    %% output to siginfo
    sig_info(k-1,:) = {gain, baseline, unit,...
        adcresolution, adczero, initialvalue,...
        checksum, blocksize, desc};
    
end

end
