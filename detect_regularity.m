function [ qrs_final, sqi_ecg, sqi_abp, ann_jqrs, ann_gqrs ] = detect_regularity(data, header, fs, opt_input)
%[ beat, sqi ] = detect_sqi(DATA, HEADER, FS) detects QRS complexes in the given
% matrix of data. HEADER must contain signal names which map to the list below.
% FS must contain a numeric sampling frequency.

% This function uses an estimate of signal quality to switch between
% signals.

%   DATA - NxD matrix of N samples with D signals. Each signal should correspond
%       to a signal name in HEADER. Do not include time as a singal.
%   HEADER - 1xD cell array of strings, containing the signal name (e.g. 'ECG', see below)
%   FS  - scalar double, containing the sampling frequency

%	This QRS detector fuses beats detected on the ECG and the ABP waveforms
%
%	HEADER	- the following are acceptable signal names
%			ECG	- electrocardiogram
%			ABP	- arterial blood pressure
%			ART	- arterial blood pressure
%			SV	- stroke volume
%			PPG	- photoplethysmogram

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
if nargin<3
    error('Function requires at least 3 inputs.');
end

% check data types
if isnumeric(data)==0
    error('Data must be an NxD matrix of N samples and D signals');
end
if isnumeric(fs)==0
    error('Third input must be a 1x1 integer representing the sampling frequency.');
end
if ischar(header) && size(data,2)==1
    header = {header}; % convert to 1x1 cell array of strings
end

% check that data/header are consistent in size
if size(data,2) ~= size(header,2)
    error('Each signal in DATA must have a corresponding element in HEADER, and vice versa.');
end

%% PRE-GAME
% if true, saves detections to WFDB format annotation files
SAVE_STUFF = 0;
recordName = ['TMP_' datestr(now,'dd-mm-yyyy-HHMMSSFFF')];

% Flag which marks suspected pacing, as indicated by large ABP delays
SUSPECTED_PACING = 0;

% Flag which turns off SV/PPG detectors
% these detectors are currently not used in the 'sqi' approach
ENABLE_OTHER_DETECTORS = 0;

[ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header);

%% Algorithm parameters which can be tuned by the user
if nargin<4 || ~isstruct(opt_input)    
    [ opt ] = setOptions;
else
    [ opt ] = setOptions(opt_input);
end
opt.LG_REC = size(data,1) ./ fs; % length of the record in seconds
opt.N_WIN = ceil(opt.LG_REC/opt.REG_WIN); % number of windows in the signal

%% Simple mode
if opt.SIMPLEMODE==1
    % we only use the first ECG and first ABP
    if numel(idxECG)>1
        idxECG = idxECG(1);
    end
    
    if numel(idxABP)>1
        idxABP = idxABP(1);
    end
end

%% Inform the user what signals will be used
fprintf('Using %d ECG and %d ABP leads. ',numel(idxECG), numel(idxABP));
if ENABLE_OTHER_DETECTORS == 0
    fprintf('Not using PPG/SV.\n');
else
    fprintf('Using PPG/SV if available.\n');
end
%% INITIALISE
% Pre-allocate space
% we will subsample cells to only contain relevant signals later

[N,M] = size(data);
ann_jqrs    = cell(1,M);
ann_gqrs    = cell(1,M);
sqi_ecg     = cell(1,M);
ppg         = cell(1,M);
sv          = cell(1,M);
sqi_bp      = cell(1,M);
abp         = cell(1,M);
sqi_abp      = cell(1,M);
abp_delay   = nan(1,M);
ppg_delay   = nan(1,M);
sv_delay    = nan(1,M);


%% Write the data out to wfdb
%FIXME: this is required because gqrs/wabp run on .dat/.hea files
% ideally, we would have these functions run on matlab data directly
wrsamp((transpose(1:N)-1),data*10000,recordName,fs,10000,'16+24');

%% ECG PEAK DETECT
if ~isempty(idxECG)
    for m=idxECG
        fprintf('Lead %d of %d. Signal %d (%s).\n',...
            find(m==idxECG,1), numel(idxECG), m, header{m});
        
        %=== call jqrs
        fprintf('\tRunning jqrs... ');
        ann_jqrs{m} = run_qrsdet_by_seg_ali(data(:,m),fs,opt);
        fprintf('done.\n');
        
        %=== call gqrs
        % usage: gqrs -r recordName -s signalIndex
        fprintf('\tRunning gqrs... ');
        gqrs(recordName,[],[],m,[],['gqrs' num2str(m)]);
        
        % load in gqrs
        try
            ann_gqrs{m} = rdann(recordName,['gqrs' num2str(m)]);
        catch
            % rdann sometimes crashes when loading empty qrs annotations
            ann_gqrs{m} = [];
        end
        fprintf('done.\n');
        
        %=== delete gqrs' output
        if SAVE_STUFF==0
            delete([recordName '.gqrs' num2str(m)]);
        end
        
        %=== convert to time
        ann_jqrs{m} = ann_jqrs{m}(:) ./ fs;
        ann_gqrs{m} = ann_gqrs{m}(:) ./ fs;
    end
end % end 'if ECG exists' segment

%% ABP (if present)
if ~isempty(idxABP)
    for m=idxABP
        switch opt.ABPMethod
            case 'delineator'
                [onsetp,peakp] = delineator(data(:,m),fs);
                abp{m} = corrDelineator(data(:,m),peakp,onsetp,fs,1);
                abp{m} = abp{m}(:); % enforce column vector
            otherwise % default is wabp for unrecognised strings
                wabp(recordName,[],[],[],m);
                try
                    abp{m} = rdann(recordName,'wabp');
                catch
                    abp{m} = [];
                end
        end
        
        % convert to seconds
        abp{m} = abp{m} / fs;
        
        switch opt.DELAYALG
            case {'crosscorr','cc'}
                %=== use cross correlation to map ABP back to ECG
                if isempty(idxECG)
                    % if there is no ECG, we can't do autocorrelation
                    % in this case, use a fixed delay of 200ms
                    abp_delay(m)=0.2;
                else
                    % we use the first ECG lead
                    % in theory we could select the highest quality lead
                    % but this is not done here
                    lagVector = mapABPtoECGcc(data(:,idxECG(1)), data(:,m), fs);
                    if numel(lagVector)>0; lagVector(isnan(lagVector)) = []; end
                    if numel(lagVector)>0; lagVector(lagVector==0) = []; end
                    if isempty(lagVector)
                        abp_delay(m) = 0.2;
                    else
                        lagVector = sort(lagVector);
                        abp_delay(m) = lagVector( ceil(numel(lagVector)/2) ) / fs; % take median
                    end
                end
                abp{m} = abp{m} - abp_delay(m);
            otherwise
                %=== we loop through the ECGs - if we find a good one, we use
                %it, otherwise, we continue to the next
                for l=idxECG
                    % this function searches for a one minute segment where
                    % the ECG and ABP annotations alternate
                    % if this segment exists, it's likely to be good
                    % quality.
                    [abp{m}, sqi_bp{m}, abp_delay(m), foundMatch] = mapABPtoRR(abp{m}, ann_gqrs{l},  sqi_bp{m});
                    
                    % this means we found a decent segment to get rrDiff
                    if foundMatch==true
                        break;
                    end
                end
        end
        
        %=== output to file
        if SAVE_STUFF == 1 && ~isempty(abp{m})
            wrann(recordName,['wabpmapped' num2str(m)],ceil(abp{m}*fs),[],[],[],[]);
            movefile([recordName '.wabp'],[recordName '.wabp' num2str(m)]);
        else
            if exist([recordName '.wabp'],'file')==2
                delete([recordName '.wabp']);
            end
        end
    end
    if opt.USE_PACING==1
        if strcmp(opt.DELAYALG,'map') && min(abp_delay) > 0.4
            SUSPECTED_PACING = 1;
            % undo the delay on ABP peaks, and replace it with default 0.2
            for m=idxABP
                abp{m} = abp{m} + abp_delay(m) - 0.2;
                abp_delay(m) = 0.2;
            end
        end
    end
end

%% SV (if present)
if ENABLE_OTHER_DETECTORS == 1 && ~isempty(idxSV)
    for m=idxSV
        sv{m} = C2014_SVDetector(data(:,m),fs);
        sv{m} = sv{m}(:) ./ fs;
        
        if ~isempty(ann_gqrs)
            [sv{m},dummy,sv_delay(m)] = mapSVtoRR(sv{m}, ann_gqrs{1});
        else
            sv_delay(m) = 0.2;
            sv{m} = ceil(sv{m} - sv_delay(m)*fs);
        end
        
        %=== output to file
        if SAVE_STUFF == 1 && ~isempty(sv{m})
            wrann(recordName,['sv' num2str(m)],sv{m},[],[],[],[]);
        end
        
    end
end

%% PPG (if present)
if ENABLE_OTHER_DETECTORS == 1 && ~isempty(idxPPG)
    for m=idxPPG
        ppg{m} = run_ppgdet_by_seg(data(:,m),fs,15,0.4,'MECG');
        ppg{m} = ppg{m}(:) ./ fs;
        ppg{m} = mapSVtoRR(ppg{m}, ann_jqrs{1});   %-- yes the interval seems to be the same as SV
        
        if ~isempty(ann_gqrs)
            [ppg{m},dummy,ppg_delay(m)] = mapSVtoRR(ppg{m}, ann_gqrs{1});
        else
            ppg_delay(m) = 0.2;
            ppg{m} = ceil(ppg{m} - ppg_delay(m)*fs);
        end
        
        %=== output to file
        if SAVE_STUFF == 1 && ~isempty(ppg{m})
            wrann(recordName,['ppg' num2str(m)],ppg{m},[],[],[],[]);
        end
    end
end

%% lead switching using regularity of the RR interval
qrs_final = cell(opt.N_WIN,1);
% put all the annotations into a single cell array
qrs_out = [ann_jqrs(idxECG),ann_gqrs(idxECG),...
    ppg(idxPPG), abp(idxABP), sv(idxSV)];
M = numel(qrs_out);

for w=1:opt.N_WIN
    curr_qrs = cell(1,M);
    ww = w*opt.REG_WIN; % in seconds
    % "SMI" is our regularity index
    % lower values indicate more regular RR interval time series
    % it's simply the standard deviation of the RR interval histogram
    % we default it to 100 - this is a *very* large value
    SMI = 100*ones(1,numel(curr_qrs)+1);
    for m=1:M
        % first, subselect a window of data
        curr_qrs{m} = qrs_out{m}(qrs_out{m}>ww-opt.REG_WIN & qrs_out{m}<=ww);
        
        % if we have more than 2 peaks, calculate std(std(RR intervals))
        if numel(curr_qrs{m})>2
            SMI(m) = assess_regularity(curr_qrs{m}-(ww-opt.REG_WIN),0,1,0.96,opt.REG_WIN,0);
        end
    end
    
    %FIXME: if all regularities are poor, we should probably output no
    %annotations rather than garbage
    [minSMI,idxMin] = min(SMI);
    qrs_final{w} = curr_qrs{idxMin};
    
    %=== we must remove double detections at the boundaries
    % this can occur if we switch from one lead to another, and when we
    % switched, two annotations occured on either side of the boundary
    % say we switched from lead 1 to lead 2 at a boundary |
    % x  x  x|  x
    %  o  o  |o  o
    % .. fused into:
    % x  x  x|o  o
    
    % the below code fixes it to:
    % x  x  x|   o
    
    % takes advantage of the heart's refractory period being 250ms
    if w>1 && ~isempty(qrs_final{w-1}) && ~isempty(qrs_final{w})
        lastQRS = qrs_final{w-1}(end);
        if abs(lastQRS - qrs_final{w}(1)) < 0.25 % within refractory
            % -> assume it is a double detection and remove
            qrs_final{w}(1) = [];
        end
    end
end

qrs_final = vertcat(qrs_final{:});
if ~isempty(qrs_final)
    qrs_final = round(qrs_final(:)*fs);
end

if ~isempty(qrs_final) && SAVE_STUFF==1
    %=== write out to file
    wrann(recordName,'qrs',qrs_final,[],[],[],[]);
end

if SAVE_STUFF==0
    % clean up folder
    delete([recordName '.dat']);
    delete([recordName '.hea']);
end
end


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


function [ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header)

bp_ind1 = cellfun(@any, regexpi(header,'bp','once'));
bp_ind2 = cellfun(@any, regexpi(header,'art','once'));
bp_ind3 = cellfun(@any, regexpi(header,'pressure','once'));

%=== retain order in ABP leads, we may prefer earlier ones for speed
bp_ind1 = find(bp_ind1==1);
bp_ind2 = find(bp_ind2==1);
bp_ind3 = find(bp_ind3==1);
bp_ind2 = setdiff(bp_ind2,bp_ind1);
bp_ind3 = setdiff(bp_ind3,[bp_ind2(:)',bp_ind1(:)']);
idxABP = [bp_ind1(:)',bp_ind2(:)',bp_ind3(:)'];

% Search for other signals
idxECG = cellfun(@any, regexpi(header,'ecg','once'));
idxECG = idxECG | cellfun(@any, regexpi(header,'ekg','once'));
idxECG = idxECG | ismember(header,...
    {'I','II','III','AVR','AVL',...
    'AVF','V','V1','V2','V3','V4',...
    'V5','V6','MCL1','MCL2','MCL3',...
    'MCL4','MCL5','MCL6','aVR','aVL','aVF'});
% lowercase letters
idxECG = idxECG | ismember(header,...
    lower({'I','II','III','AVR','AVL',...
    'AVF','V','V1','V2','V3','V4',...
    'V5','V6','MCL1','MCL2','MCL3',...
    'MCL4','MCL5','MCL6','aVR','aVL','aVF'}));
% uppercase letters
idxECG = idxECG | ismember(header,...
    upper({'I','II','III','AVR','AVL',...
    'AVF','V','V1','V2','V3','V4',...
    'V5','V6','MCL1','MCL2','MCL3',...
    'MCL4','MCL5','MCL6','aVR','aVL','aVF'}));

idxECG = find(idxECG);

idxSV = cellfun(@any, regexpi(header,'sv','once'));
idxSV = find(idxSV==1);

idxPPG = cellfun(@any, regexpi(header,'ppg','once'));
idxPPG = idxPPG | cellfun(@any, regexpi(header,'pleth','once'));
idxPPG = find(idxPPG==1);

% ensure all vectors are row vectors
% this facilitates their use in a for loop, i.e. for m=idxECG
idxECG = idxECG(:)';
idxABP = idxABP(:)';
idxPPG = idxPPG(:)';
idxSV = idxSV(:)';
end


function [ opt ] = setOptions(varargin)
%=== parameters that define the window for the bSQI check on the ECG
opt_default.SIZE_WIND = 10;
opt_default.HALF_WIND = opt_default.SIZE_WIND/2;

% take the median SQI using X nearby values
% this is used only for the ECG SQI
opt_default.LG_MED = 3;
% so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows

% how frequently to check the SQI for switching
% i.e., if REG_WIN = 1, then we check the signals every second to switch
opt_default.REG_WIN = 1;

% the width, in seconds, used when comparing peaks in the F1 based ECG SQI
opt_default.THR = 0.150;

% the SQI threshold
%   if the lead SQI is higher than this we use this signal
%   if the lead SQI is lower than this, we use the next signal
%   this an ordered comparison, so we usually default to the ECG signal
%   (the first column/signal present in the data)
opt_default.SQI_THR = 0.8;
opt_default.USE_PACING = 1; % flag turning on/off the pacing detection/correction

% ABP peak detection method
%   options:
%       wabp
%       delineator
opt_default.ABPMethod = 'wabp';

opt_default.SIMPLEMODE = 0; % only use the first of ABP/ECG

%=== jqrs parameters
opt_default.JQRS_THRESH = 0.3;
opt_default.JQRS_REFRAC = 0.25;
opt_default.JQRS_INTWIN_SZ = 7;
opt_default.JQRS_WINDOW = 15;

opt_default.DELAYALG = 'map';

if nargin==0
    opt = opt_default;
    return;
elseif nargin==1
    opt = varargin{1};
else
    error('Incorrect number of inputs.');
end

if isfield(opt,'DELAYALG')
    if ischar(opt.DELAYALG)~=1 || any(ismember(opt.DELAYALG,{'map','crosscorr','cc'}))==0
        fprintf('Delay algorithm name unrecognised - using default peak mapping.\n');
        opt.DELAYALG = 'map';
    end
end

% if input options are given, we update the default opt with them
if nargin>0 && isstruct(opt)
    fn = fieldnames(opt);
    fn_default = fieldnames(opt_default);
    for f=1:numel(fn)
        if ismember(fn{f},fn_default)
            opt_default.(fn{f}) = opt.(fn{f});
        end
    end
end
opt = opt_default;

end
