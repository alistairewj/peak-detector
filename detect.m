function [ qrs_final, sqi_ecg, ann_jqrs, ann_gqrs ] = detect(recordName, FUSEALG, DELAYALG, opt_input)
%[ beat, sqi ] = detect(RECORDNAME) detects QRS complexes in the given
%record. RECORDNAME must be in a WFDB compatible file format.

%	This QRS detector fuses beats detected on the ECG and the ABP waveforms
%
%	signals	- the following are acceptable signal names
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

%% PRE-GAME
% if true, saves detections to annotation files and re-saves the data to a mat file
SAVE_STUFF = 0;
if nargin<2
FUSEALG = 'sqi';
end
if nargin<3
DELAYALG = 'map';
end

if ischar(FUSEALG)~=1 || any(ismember(FUSEALG,{'sqi','regularity'}))==0
    fprintf('Fusion algorithm name unrecognised - using default SQI switching.\n');
    FUSEALG = 'sqi';
end

if ischar(DELAYALG)~=1 || any(ismember(DELAYALG,{'map','crosscorr','cc'}))==0
    fprintf('Delay algorithm name unrecognised - using default peak mapping.\n');
    DELAYALG = 'map';
end


% Flag which marks suspected pacing, as indicated by large ABP delays
SUSPECTED_PACING = 0;

% Flag which turns on/off SV/PPG detectors
% these detectors are currently not used in the 'sqi' approach
switch FUSEALG
    case 'regularity'
        ENABLE_OTHER_DETECTORS = 1;
    otherwise
        ENABLE_OTHER_DETECTORS = 0;
end

%% LOAD DATA
[ data, header, fs ] = loadData(recordName);
[ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header);

%% ECG sqi parameters
opt = struct();
if nargin<4 || ~isstruct(opt_input)    
    %=== parameters that define the window for the bSQI check on the ECG
    opt.SIZE_WIND = 10;
    opt.HALF_WIND = opt.SIZE_WIND/2;
    
    % take the median SQI using X nearby values
    % this is used only for the ECG SQI
    opt.LG_MED = 3;
    % so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows
    
    % how frequently to check the SQI for switching
    % i.e., if REG_WIN = 1, then we check the signals every second to switch
    opt.REG_WIN = 1;
    
    % the width, in seconds, used when comparing peaks in the F1 based ECG SQI
    opt.THR = 0.150;
    
    % the SQI threshold
    %   if the lead SQI is higher than this we use this signal
    %   if the lead SQI is lower than this, we use the next signal
    %   this an ordered comparison, so we usually default to the ECG signal
    %   (the first column/signal present in the data)
    opt.SQI_THR = 0.8;
    opt.USE_PACING = 1; % flag turning on/off the pacing fix
    
    % ABP peak detection method
    %   options:
    %       wabp
    %       delineator
    opt.ABPMethod = 'wabp';
    
    opt.SIMPLEMODE = 0; % only use the first of ABP/ECG
    
    %=== jqrs parameters
    opt.JQRS_THRESH = 0.3;
    opt.JQRS_REFRAC = 0.25;
    opt.JQRS_INTWIN_SZ = 7;
    opt.JQRS_WINDOW = 15;
else
    [ opt ] = setOptions(opt_input);
end
opt.LG_REC = size(data,1) ./ fs; % length of the record in seconds
opt.N_WIN = ceil(opt.LG_REC/opt.REG_WIN); % number of windows in the signal
LINUX_FLAG = isunix;

opt.SIMPLEMODE=1;

if opt.SIMPLEMODE==1
    % we only use the first ECG and first ABP
    if numel(idxECG)>1
        idxECG = idxECG(1);
    end
    
    if numel(idxABP)>1
        idxABP = idxABP(1);
    end
end

fprintf('Using %d ECG and %d ABP leads.\n',numel(idxECG), numel(idxABP));
if ENABLE_OTHER_DETECTORS == 0
    fprintf('Not using PPG/SV\n');
else
    fprintf('Using PPG/SV\n');
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
abpsqi      = cell(1,M);
abp_delay   = nan(1,M);
ppg_delay   = nan(1,M);
sv_delay    = nan(1,M);

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
        if LINUX_FLAG
            system(['gqrs -r ' recordName ' -f 0 -o qrs -s ' num2str(m-1)]);
        else
            gqrs(recordName,[],[],m,[],'qrs');
        end
        
        % load in gqrs
        if LINUX_FLAG == 1
        system(['rdann -r ' recordName ' -a qrs -x | awk ''{print $1}'' > tmpout']);
        fp = fopen('tmpout','r');
        tmp = textscan(fp,'%s','delimiter','\n');
        tmp = tmp{1};
        tmp = str2double(tmp);
        ann_gqrs{m} = round(tmp * fs)+1; % convert 0-1 indexing
        else
            
        try
            ann_gqrs{m} = rdann(recordName,'qrs');
        catch
            % rdann sometimes crashes when loading empty qrs annotations
            ann_gqrs{m} = [];
        end
        end
        fprintf('done.\n');
        
        %=== output to file
        if SAVE_STUFF == 1
            movefile([recordName '.qrs'],[recordName '.gqrs' num2str(m)]);
        else
            delete([recordName '.qrs']);
        end
        
        %=== convert to time
        ann_jqrs{m} = ann_jqrs{m}(:) ./ fs;
        ann_gqrs{m} = ann_gqrs{m}(:) ./ fs;
        
    end
end % end 'if ECG exists' segment
% ann_gqrs    = ann_gqrs(idxECG);
% ann_jqrs    = ann_jqrs(idxECG);
% sqi_ecg     = sqi_ecg(idxECG);


%% ABP (if present)
if ~isempty(idxABP)
    for m=idxABP
        switch opt.ABPMethod
            case 'delineator'
                [onsetp,peakp] = delineator(data(:,m),fs);
                abp{m} = corrDelineator(data(:,m),peakp,onsetp,fs,1);
                abp{m} = abp{m}(:); % enforce column vector
            otherwise % default is wabp for unrecognised strings
                if LINUX_FLAG
                    system(['wabp -r ' recordName ' -s ' num2str(m-1)]);
                else
                    wabp(recordName,[],[],[],m);
                end
                
                if LINUX_FLAG == 1
                system(['rdann -r ' recordName ' -a wabp -x | awk ''{print $1}'' > tmpout']);
                fp = fopen('tmpout','r');
                tmp = textscan(fp,'%s','delimiter','\n');
                tmp = tmp{1};
                tmp = str2double(tmp);
                abp{m} = round(tmp * fs)+1; % convert 0-1 indexing
                else
                    try
                    abp{m} = rdann(recordName,'wabp');
                    catch
                        abp{m} = [];
                    end
                end
%                 try
%                     abp{m} = rdann(recordName,'wabp');
%                 catch
%                     % rdann sometimes crashes with empty annotations
%                     abp{m} = [];
%                 end
        end
        
        
        if strcmp(FUSEALG,'sqi')
            %=== get SQI for ABP signals
            [ sqi_bp{m}, header_sq, sqi_bp_values, header_sqi ] = calcABPSQI(data(:,m), abp{m}, fs);
        end
        abp{m} = abp{m} / fs;
        
        switch DELAYALG
            case {'crosscorr','cc'}
                %=== use cross correlation to map ABP back to ECG
                if isempty(idxECG)
                    abp_delay(m)=0.2; % hardcoded delay if no ECG is available
                else
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
    if strcmp(DELAYALG,'map') && min(abp_delay) > 0.4 
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

%% switching
qrs_final = cell(opt.N_WIN,1);
switch FUSEALG
    case 'regularity'
        %% combine all the annotations
        qrs_out = [ann_jqrs(idxECG),ann_gqrs(idxECG),...
            ppg(idxPPG), abp(idxABP), sv(idxSV)];
        M = numel(qrs_out);
        % regularity!!
        for w=1:opt.N_WIN
            curr_qrs = cell(1,M);
            ww = w*opt.REG_WIN; % in seconds
            SMI = 100*ones(1,numel(curr_qrs)+1);
            for m=1:M
                curr_qrs{m} = qrs_out{m}(qrs_out{m}>ww-opt.REG_WIN & qrs_out{m}<=ww);
                
                if numel(curr_qrs{m})>2
                    SMI(m) = assess_regularity(curr_qrs{m}-(ww-opt.REG_WIN),0,1,0.96,opt.REG_WIN,0);
                end
            end
            
            [minSMI,idxMin] = min(SMI);
            qrs_final{w} = curr_qrs{idxMin};
            
            %=== we must remove double detections at the boundaries
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
        
    otherwise % SQI switching
        %% ECG LEAD WISE SQI
        for m=idxECG
            fprintf('\tCalculating windowed sqi... ');
            if ~isempty(ann_gqrs{m}) && ~isempty(ann_jqrs{m})
                %=== calculate the SQI for each ECG lead
                % all tsqi (start time of SQI window) will be identical, so
                % no need for a cell to store these
                [ sqi_ecg{m}, tsqi ] = ecgsqi( ann_gqrs{m}, ann_jqrs{m}, opt );
            else
                sqi_ecg{m} = zeros(opt.N_WIN,1);
                tsqi = 0:opt.REG_WIN:opt.LG_REC;
                tsqi(tsqi==opt.LG_REC) = [];
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
                abpsqi{m} = mean(sqi_bp{m}(:,1));
            %=== we are switching less frequently than our SQI win size
            % so we can calculate window by window SQIs for ABP easily
            elseif xi(3) >= xi(2)
                xi = xi(:);
                idxLast = find(xi>opt.LG_REC,1);
                if mod(idxLast,2)==0; xi = xi(1:idxLast); end
                %=== first we create an index, idxMap
                % this represents which window each beat belongs to
                % we can then use idxMap to determine if a beat is in window w
                [tmp,idxMap] = histc(abp{m},xi);
                tmp = tmp(1:max(idxMap));
                tmp2 = accumarray(idxMap(idxMap~=0),sqi_bp{m}(idxMap~=0,1)) ./ tmp;
                idxKeep = unique(idxMap);
                idxKeep(idxKeep==0) = [];
                abpsqi{m} = zeros(size(xi));
                abpsqi{m}(idxKeep) = tmp2(idxKeep);
                abpsqi{m} = abpsqi{m}(1:2:end);
            else
                %=== we have more than 1 switch event per SQI window
                % this prevents directy application of histc, have to be a
                % bit more clever
                
                % create windows for the SQI
                tabpsqi = bsxfun(@plus, ...
                    0:opt.SIZE_WIND:opt.LG_REC+opt.SIZE_WIND-0.01, ...
                    transpose(0:opt.REG_WIN:opt.SIZE_WIND-0.01)); 
                % -0.01 prevents duplicate window at bottom of xi_all
                
                % each column in xi_all moves through the signal
                % each row is a different starting point from 0:SIZE_WIN
                % so each row allows us to calculate the SQI using histc
                % then we move to the next row to do the same
                % then we vertically concatenate it all together
                
                abpsqi{m} = zeros(opt.SIZE_WIND/opt.REG_WIN,size(tabpsqi,2)-1);
                for d=1:size(tabpsqi,1)
                    xi = tabpsqi(d,:);
                    % idxMap finds which window each abp beat corresponds
                    [tmp,idxMap] = histc(abp{m},xi);
                    % accumarray sums the ABP SQIs with the same idxMap
                    % dividing by tmp makes abpsqi{m} the mean SQI
                    abpsqi{m}(d,1:max(idxMap)) = accumarray(idxMap(idxMap~=0),sqi_bp{m}(idxMap~=0,1)) ./ tmp(1:max(idxMap));
                end
                abpsqi{m} = abpsqi{m}(:);
                
                %=== values of 'NaN' are due to a division by 0
                % these occur if there are no beats in the given window
                % we should probably call these segments bad quality!
                abpsqi{m}(isnan(abpsqi{m})) = 0;
                
                %=== now we create the final abp times
                tabpsqi = tabpsqi(:,1:end-1); % histc always excludes last col
                tabpsqi = tabpsqi(:);
                
                %=== now we delete the indices after the end of the rec
                abpsqi{m}(tabpsqi>=opt.LG_REC) = [];
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
        qrs_in = [ann_gqrs(idxECG), abp(idxABP)];
        sqi_in = [sqi_ecg(idxECG), abpsqi(idxABP)];
        
        % use a delay of 100ms for the ECG, as it is used in the switching
        % function to check for missed beats on SQI switching transitions
        abp_delay = [0.1*ones(1,numel(idxECG)),abp_delay(idxABP)];
        %=== perform switching
        if numel(qrs_in)==1
            qrs_final = round(qrs_in{:});
        else
            for m=1:numel(sqi_in)
                if numel(sqi_in{m})<opt.N_WIN
                    sqi_in{m} = [sqi_in{m};zeros(opt.N_WIN - numel(sqi_in{m}),1)];
                elseif numel(sqi_in{m})>opt.N_WIN
                    sqi_in{m} = sqi_in{m}(1:opt.N_WIN);
                end
            end
            qrs_final = sqi_switching(qrs_in,sqi_in,tsqi,opt.SQI_THR,abp_delay,0);
        end
        
        %=== write out annotation with a hack for rounding
	qrs_final = qrs_final(:)*fs;
	qrs_final = round(qrs_final*10);
        qrs_final = round(qrs_final/10);
end

if ~isempty(qrs_final)
    %=== write out to file
    wrann(recordName,'qrs',qrs_final,[],[],[],[]);
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


function [ opt ] = setOptions(opt)

opt_default.SQI_THR = 0.8; % default threshold for switching
opt_default.THR = 0.15; % window for matching two peaks
opt_default.LG_MED = 3; % take the median SQI across X seconds
opt_default.SIZE_WIND = 10;
opt_default.REG_WIN = 1; % one window per second
opt_default.HALF_WIND = opt_default.SIZE_WIND/2;
opt_default.LG_REC = NaN;  % length of the record in seconds
opt_default.N_WIN = NaN; % number of windows in the signal
opt_default.USE_PACING = 1; % use pacing detection and correction
opt_default.ABPMethod = 'wabp';
opt_default.SIMPLEMODE = 0;

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
