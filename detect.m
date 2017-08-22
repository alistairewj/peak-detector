function [ qrs, sqi ] = detect(recordName, alg, opt)
%[ qrs, sqi ] = detect(RECORDNAME) detects QRS complexes in the
%WFDB-readable record RECORDNAME. The QRS detections and signal quality
%estimates are ouput.

% Dependencies:
%
%       1) This function requires the WFDB Toolbox for MATLAB and Octave.
%          For information on how to install the toolbox, please see:
%              http://www.physionet.org/physiotools/matlab/wfdb-app-matlab/
%
%       2) The Toolbox is supported only on 64-bit MATLAB 2013a and 2013b
%          and on GNU Octave 3.6.4 and later, on Linux, Mac OS X, and
%          Windows.

%% CHECK INPUTS ARE VALID

if nargin<1
    error('Function requires at least 1 inputs.');
end

% check for the data/header files
if exist([recordName '.hea'],'file') ~= 2
    error('detect:noHeaderFile',...
        ['Unable to locate header file ''' recordName '''.hea'])
end
if exist([recordName '.dat'],'file') ~= 2 && exist([recordName '.mat'],'file') ~= 2 
    error('detect:noDataFile',...
        ['Unable to locate data file ''' recordName '''.dat'])
end

if nargin<2
    alg = 'sqi';
elseif ~ischar(alg)
    error('detect:incorrectAlgorithm',...
        'Switching algorithm should be a string: either ''sqi'' or ''regularity''.')
elseif ~strcmpi(alg,'sqi') && ~strcmpi(alg,'regularity')
    error('detect:incorrectAlgorithm',...
        'Switching algorithm should be ''sqi'' or ''regularity''')
end

% if no options input, create empty structure
% option settings are checked in the detect subfunctions
if nargin<3
    opt = struct();
end

%% LOAD DATA
[data,fs] = rdsamp(recordName);
[siginfo,fs] = wfdbdesc(recordName);

% extract info from structure output by wfdbdesc
header = arrayfun(@(x) x.Description, siginfo, 'UniformOutput', false);

%% run either SQI algorithm or regularity algorithm
if strcmpi(alg,'regularity')
    % run heart rate regularity based switching
    [ qrs, sqi, qrs_comp, qrs_header ] = detect_regularity(data, header, fs, opt);
else
    % run SQI based switching (it's the only other option)
    [ qrs, sqi, qrs_comp, qrs_header ] = detect_sqi(data, header, fs, opt);
end

end

% optionally we can load the data using wfdb2mat, which is much faster
function [ data, header, fs ] = loadData(recordName)
[d,config] = wfdbloadlib; % set wfdb library

if isunix
    [err,res] = system(['wfdb2mat -r ' recordName ' -H']);
else
    err=1;
end
if err~=0 % if above call fails, use rdsamp
    %=== only load data using rdsamp .. and hope it's in physical units
    [~,data,fs] = rdsamp(recordName);
    
    %=== get bp signal index
    siginfo = wfdbdesc(recordName);
    header = cell(1,numel(siginfo));
    for k=1:numel(siginfo)
        header{k} = siginfo(k).Description;
    end    
else
    % parse the output text for signal FS and gain
    res = regexp(res,'\n','split');
    fs = find(strncmp(res,'Sampling frequency',18)==1,1);
    fs = res{fs}(20:regexp(res{fs},'Hz','once')-1);
    fs = str2double(strtrim(fs));
    
    % isolate the rows with signals
    res = res(find(strncmp(res,'Row',3)==1,1)+1:end-3);
    if strcmp(res(end),'') == 1; res = res(1:end-1); end
    N_SIG = numel(res);
    header = repmat({'NULL'},N_SIG,1);
    gain = ones(1,N_SIG);
    base = zeros(1,N_SIG);
    
    for m=1:N_SIG
        tmp = res{m};
        tmp = regexp(tmp,'\t','split');
        if numel(tmp)>3
            header{m} = tmp{2};
            gain(m) = str2double(tmp{3});
            base(m) = str2double(tmp{4});
        end
    end
    
    data = load([recordName 'm']);
    data = transpose(data.val);
    
    %=== apply gain // offset
    for m=1:N_SIG
        data(:,m) = (data(:,m) - base(m))/gain(m);
    end
end


end
