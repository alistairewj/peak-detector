function [ qrs, smi, qrs_comp, qrs_header ] = detect_regularity(data, header, fs, opt)
%[ QRS ] = detect_sqi(DATA, HEADER, FS) detects QRS complexes in the given
% matrix of data. HEADER must contain signal names which map to the list below.
% FS must contain a numeric sampling frequency.
%
%[ QRS, QRS_COMP, QRS_HEADER ] = detect_sqi(DATA, HEADER, FS) also returns
%the constitutent QRS detectors (QRS_COMP) and a header describing each
%detector (QRS_HEADER).
%
% This function uses an estimate of signal quality to switch between
% signals.
%
%   DATA - NxD matrix of N samples with D signals. Each signal should correspond
%       to a signal name in HEADER. Do not include time as a singal.
%   HEADER - 1xD cell array of strings, containing the signal name (e.g. 'ECG', see below)
%   FS  - scalar double, containing the sampling frequency
%
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

%% option setting

% check if options are input
if nargin<4
    % default options
    [ opt ] = setDetectOptions;
else
    if ~isstruct(opt)
        error('detect:invalidOptions',...
            'Second argument should be a structure.');
    else
        [ opt ] = setDetectOptions(opt);
    end
end

opt.LG_REC = size(data,1) ./ fs; % length of the record in seconds
opt.N_WIN = ceil(opt.LG_REC/opt.REG_WIN); % number of windows in the signal


%% PRE-GAME
% if true, saves detections to WFDB format annotation files
SAVE_STUFF = opt.SAVE_STUFF;
recordName = ['TMP_' datestr(now,'dd-mm-yyyy-HHMMSSFFF')];

[ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header);

if numel(fs) == 1
    % assume all signals have the same sampling frequency
    fs = repmat(fs,1,numel(header));
end

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

% Flag which turns off SV/PPG detectors
% these detectors are currently not used in the 'sqi' approach
fprintf('Using %d ECG and %d ABP leads. ',numel(idxECG), numel(idxABP));
if opt.ENABLE_OTHER_DETECTORS == 0
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
ppg         = cell(1,M);
sv          = cell(1,M);
abp         = cell(1,M);
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
        ann_jqrs{m} = run_qrsdet_by_seg_ali(data(:,m),fs(m),opt);
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
        ann_jqrs{m} = ann_jqrs{m}(:) ./ fs(m);
        ann_gqrs{m} = ann_gqrs{m}(:) ./ fs(m);
    end
end % end 'if ECG exists' segment

%% ABP (if present)
if ~isempty(idxABP)
    for m=idxABP
        switch opt.ABPMethod
            case 'delineator'
                [onsetp,peakp] = delineator(data(:,m),fs(m));
                abp{m} = corrDelineator(data(:,m),peakp,onsetp,fs(m),1);
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
        abp{m} = abp{m} / fs(m);
        
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
                    lagVector = mapABPtoECGcc(data(:,idxECG(1)), data(:,m), fs(m));
                    if numel(lagVector)>0; lagVector(isnan(lagVector)) = []; end
                    if numel(lagVector)>0; lagVector(lagVector==0) = []; end
                    if isempty(lagVector)
                        abp_delay(m) = 0.2;
                    else
                        lagVector = sort(lagVector);
                        abp_delay(m) = lagVector( ceil(numel(lagVector)/2) ) / fs(m); % take median
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
            wrann(recordName,['wabpmapped' num2str(m)],ceil(abp{m}*fs(m)),[],[],[],[]);
            movefile([recordName '.wabp'],[recordName '.wabp' num2str(m)]);
        else
            if exist([recordName '.wabp'],'file')==2
                delete([recordName '.wabp']);
            end
        end
    end
end

%% SV (if present)
if opt.ENABLE_OTHER_DETECTORS == 1 && ~isempty(idxSV)
    for m=idxSV
        sv{m} = C2014_SVDetector(data(:,m),fs(m));
        sv{m} = sv{m}(:) ./ fs(m);
        
        if ~isempty(ann_gqrs)
            [sv{m},dummy,sv_delay(m)] = mapSVtoRR(sv{m}, ann_gqrs{1});
        else
            sv_delay(m) = 0.2;
            sv{m} = ceil(sv{m} - sv_delay(m)*fs(m));
        end
        
        %=== output to file
        if SAVE_STUFF == 1 && ~isempty(sv{m})
            wrann(recordName,['sv' num2str(m)],sv{m},[],[],[],[]);
        end
        
    end
end

%% PPG (if present)
if opt.ENABLE_OTHER_DETECTORS == 1 && ~isempty(idxPPG)
    for m=idxPPG
        ppg{m} = run_ppgdet_by_seg(data(:,m),fs(m),15,0.4,'MECG');
        ppg{m} = ppg{m}(:) ./ fs(m);
        ppg{m} = mapSVtoRR(ppg{m}, ann_jqrs{1});   %-- yes the interval seems to be the same as SV
        
        if ~isempty(ann_gqrs)
            [ppg{m},dummy,ppg_delay(m)] = mapSVtoRR(ppg{m}, ann_gqrs{1});
        else
            ppg_delay(m) = 0.2;
            ppg{m} = ceil(ppg{m} - ppg_delay(m)*fs(m));
        end
        
        %=== output to file
        if SAVE_STUFF == 1 && ~isempty(ppg{m})
            wrann(recordName,['ppg' num2str(m)],ppg{m},[],[],[],[]);
        end
    end
end

%% lead switching using regularity of the RR interval
N_WIN_MAX = max(opt.N_WIN);

qrs = cell(N_WIN_MAX,1);
% put all the annotations into a single cell array
qrs_comp = [ann_jqrs(idxECG),ann_gqrs(idxECG),...
    ppg(idxPPG), abp(idxABP), sv(idxSV)];

%=== create the header
qrs_header = [strcat(repmat({'jqrs'}, 1, sum(~isempty(ann_jqrs(idxECG)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxECG))), 'UniformOutput', false)),...
    strcat(repmat({'gqrs'}, 1, sum(~isempty(ann_jqrs(idxABP)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxABP))), 'UniformOutput', false)),...
    strcat(repmat({'ppg'}, 1, sum(~isempty(ann_jqrs(idxPPG)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxPPG))), 'UniformOutput', false)),...
    strcat(repmat({'abp'}, 1, sum(~isempty(ann_jqrs(idxABP)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxABP))), 'UniformOutput', false)),...
    strcat(repmat({'sv'}, 1, sum(~isempty(ann_jqrs(idxSV)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxSV))), 'UniformOutput', false))];

M = numel(qrs_comp);

% "SMI" is our regularity index
% lower values indicate more regular RR interval time series
% it's simply the standard deviation of the RR interval histogram
% we default it to 100 - this is a *very* large value
smi = 100*ones(N_WIN_MAX,M+1);

for w=1:N_WIN_MAX
    curr_qrs = cell(1,M);
    ww = w*opt.REG_WIN; % in seconds
    for m=1:M
        % first, subselect a window of data
        curr_qrs{m} = qrs_comp{m}(qrs_comp{m}>ww-opt.REG_WIN & qrs_comp{m}<=ww);
        
        % if we have more than 2 peaks, calculate std(std(RR intervals))
        if numel(curr_qrs{m})>2
            smi(w,m) = assess_regularity(curr_qrs{m}-(ww-opt.REG_WIN),0,1,0.96,opt.REG_WIN,0);
        end
    end
    
    %FIXME: if all regularities are poor, we should probably output no
    %annotations rather than garbage
    [minSMI,idxMin] = min(smi(w,:));
    qrs{w} = curr_qrs{idxMin};
    
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
    if w>1 && ~isempty(qrs{w-1}) && ~isempty(qrs{w})
        lastQRS = qrs{w-1}(end);
        if abs(lastQRS - qrs{w}(1)) < 0.25 % within refractory
            % -> assume it is a double detection and remove
            qrs{w}(1) = [];
        end
    end
end

smi = smi(:,1:end-1);
qrs = vertcat(qrs{:});

if SAVE_STUFF==0
    % clean up folder
    delete([recordName '.dat']);
    delete([recordName '.hea']);
end
end
