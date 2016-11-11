function [ qrs, sqi, qrs_comp, qrs_header ] = detect_sqi(data, header, fs, opt)
%[ beat, sqi ] = detect_sqi(DATA, HEADER, FS) detects QRS complexes in the given
% matrix of data. HEADER must contain signal names which map to the list below.
% FS must contain a numeric sampling frequency.
%
%[ QRS, QRS_COMP, QRS_HEADER ] = detect_sqi(DATA, HEADER, FS) also returns
%the constitutent QRS detectors (QRS_COMP) and a header describing each
%detector (QRS_HEADER).

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

% note that we need to write out files in order to call gqrs
% this is really slow.. could do with improvement..
recordName = ['TMP_' datestr(now,'dd-mm-yyyy-HHMMSSFFF')];

% Flag which marks suspected pacing, as indicated by large ABP delays
SUSPECTED_PACING = 0;

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
%TODO: this is required because gqrs/wabp run on .dat/.hea files
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
        
        %=== get SQI for ABP signals
        [ sqi_bp{m}, header_sq, sqi_bp_values, header_sqi ] = calcABPSQI(data(:,m), abp{m}, fs(m));
        abp{m} = abp{m} / fs(m);
        
        switch opt.DELAYALG
            case {'crosscorr','cc'}
                %=== use cross correlation to map ABP back to ECG
                if isempty(idxECG)
                    abp_delay(m)=0.2; % hardcoded delay if no ECG is available
                else
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

%% ECG LEAD WISE SQI
for m=idxECG
    fprintf('\tCalculating windowed sqi on signal %d... ', m);
    if ~isempty(ann_gqrs{m}) && ~isempty(ann_jqrs{m})
        %=== calculate the SQI for each ECG lead
        % all tsqi (start time of SQI window) will be identical, so
        % no need for a cell to store these
        [ sqi_ecg{m}, tsqi ] = ecgsqi( ann_gqrs{m}, ann_jqrs{m},...
            opt.THR, opt.SIZE_WIND, opt.REG_WIN, opt.LG_MED,...
            opt.LG_REC(m), opt.N_WIN(m));
    else
        sqi_ecg{m} = zeros(opt.N_WIN(m),1);
        tsqi = 0:opt.REG_WIN:opt.LG_REC(m);
        tsqi(tsqi==opt.LG_REC(m)) = [];
        tsqi = tsqi(:);
    end
    fprintf('done.\n');
end

for m=idxABP
    %=== convert SQI from beat wise to window by window
    xi = [0:opt.REG_WIN:opt.LG_REC+(opt.REG_WIN/2);
        (0:opt.REG_WIN:opt.LG_REC+(opt.REG_WIN/2)) + opt.SIZE_WIND];
    if numel(xi)<=2
        %=== we never switch because the record is too short
        sqi_abp{m} = mean(sqi_bp{m}(:,1));
        %=== we are switching less frequently than our SQI win size
        % so we can calculate window by window SQIs for ABP easily
    elseif xi(3) >= xi(2)
        xi = xi(:);
        idxLast = find(xi>opt.LG_REC(m),1);
        if mod(idxLast,2)==0; xi = xi(1:idxLast); end
        %=== first we create an index, idxMap
        % this represents which window each beat belongs to
        % we can then use idxMap to determine if a beat is in window w
        [tmp,idxMap] = histc(abp{m},xi);
        tmp = tmp(1:max(idxMap));
        
        % if tmp has 0s, we set them to 1
        % the ABP will consequently be 0 / 1 == 0
        idxBad = tmp==0;
        if any(idxBad)
            tmp(idxBad) = 1;
        end
        
        tmp2 = accumarray(idxMap(idxMap~=0),sqi_bp{m}(idxMap~=0,1)) ./ tmp;
        idxKeep = unique(idxMap);
        idxKeep(idxKeep==0) = [];
        sqi_abp{m} = zeros(size(xi));
        sqi_abp{m}(idxKeep) = tmp2(idxKeep);
        sqi_abp{m} = sqi_abp{m}(1:2:end);
    else
        %=== we have more than 1 switch event per SQI window
        % this prevents directy application of histc, have to be a
        % bit more clever
        
        % create windows for the SQI
        tabpsqi = bsxfun(@plus, ...
            0:opt.SIZE_WIND:opt.LG_REC(m)+opt.SIZE_WIND-0.01, ...
            transpose(0:opt.REG_WIN:opt.SIZE_WIND-0.01));
        % -0.01 prevents duplicate window at bottom of xi_all
        
        % each column in xi_all moves through the signal
        % each row is a different starting point from 0:SIZE_WIN
        % so each row allows us to calculate the SQI using histc
        % then we move to the next row to do the same
        % then we vertically concatenate it all together
        
        sqi_abp{m} = zeros(opt.SIZE_WIND/opt.REG_WIN,size(tabpsqi,2)-1);
        for d=1:size(tabpsqi,1)
            xi = tabpsqi(d,:);
            % idxMap finds which window each abp beat corresponds
            [tmp,idxMap] = histc(abp{m},xi);
            
            % if tmp has 0s, we set them to 1
            % the ABP will consequently be 0 / 1 == 0
            idxBad = tmp==0;
            if any(idxBad)
                tmp(idxBad) = 1;
            end
            
            % accumarray sums the ABP SQIs with the same idxMap
            % dividing by tmp makes abpsqi{m} the mean SQI
            sqi_abp{m}(d,1:max(idxMap)) = accumarray(idxMap(idxMap~=0),sqi_bp{m}(idxMap~=0,1)) ./ tmp(1:max(idxMap));
        end
        sqi_abp{m} = sqi_abp{m}(:);
        
        %=== values of 'NaN' are due to a division by 0
        % these occur if there are no beats in the given window
        % we should probably call these segments bad quality!
        sqi_abp{m}(isnan(sqi_abp{m})) = 0;
        
        %=== now we create the final abp times
        tabpsqi = tabpsqi(:,1:end-1); % histc always excludes last col
        tabpsqi = tabpsqi(:);
        
        %=== now we delete the indices after the end of the rec
        sqi_abp{m}(tabpsqi>=opt.LG_REC(m)) = [];
    end
end

%JOs update ecg_sqi for suspected pacing beats...
% force usage of abp if good quality
if opt.USE_PACING == 1 && SUSPECTED_PACING==1
    for l=idxECG
        sqi_ecg{l} = zeros(size(sqi_ecg{l}));
    end
end

%=== set up for switching
qrs_comp = [ann_gqrs(idxECG), abp(idxABP)];
sqi = [sqi_ecg(idxECG), sqi_abp(idxABP)];


%=== create the header
qrs_header = [strcat(repmat({'gqrs'}, 1, sum(~isempty(ann_jqrs(idxECG)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxECG))), 'UniformOutput', false)),...
    %strcat(repmat({'gqrs'}, 1, sum(~isempty(ann_jqrs(idxABP)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxABP))), 'UniformOutput', false)),...
    strcat(repmat({'abp'}, 1, sum(~isempty(ann_jqrs(idxABP)))), arrayfun(@num2str, 1:sum(~isempty(ann_jqrs(idxABP))), 'UniformOutput', false))];


% use a delay of 100ms for the ECG, as it is used in the switching
% function to check for missed beats on SQI switching transitions
abp_delay = [0.1*ones(1,numel(idxECG)),abp_delay(idxABP)];

N_WIN_MAX = max(opt.N_WIN);

%=== perform switching
if numel(qrs_comp)==1
    % only one set of QRS detections is present - so no switching is
    % possible
    qrs = qrs_comp{1};
else
    % for each SQI provided
    for m=1:numel(sqi)
        % if we have fewer windows in this lead, extend it to the signal
        % size
        if numel(sqi{m})<N_WIN_MAX
            sqi{m} = [sqi{m};zeros(N_WIN_MAX - numel(sqi{m}),1)];
        elseif numel(sqi{m})>N_WIN_MAX
            % this if seems superfluous - we should never have more windows
            % than the max
            sqi{m} = sqi{m}(1:opt.N_WIN_MAX);
        end
    end
    qrs = sqi_switching(qrs_comp,sqi,tsqi,opt.SQI_THR,abp_delay,0);
end

if SAVE_STUFF==0
    % clean up folder
    delete([recordName '.dat']);
    delete([recordName '.hea']);
end
end

