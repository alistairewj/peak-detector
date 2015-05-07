function [qrs_final] = hailmary(data,header,fs,inputqrs,opt_sqi,fuse)
% Wavelet based feature signals
%
% This function generates a couple more feature signals for QRS detection.
% It relies on the ECG signals detections.
%
% Input:
%   data           Raw data containig all signals [MxN matrix]
%   header         Header file for channel info
%   fs             Sampling frequency [Hz]
%   inputqrs       Cell containing available qrs detections, so that:
%                       - inputqrs{1} = ECG qrs detections
%                       - inputqrs{2} = PPG mapped detections
%                       - inputqrs{3} = ABP mapped detections
%                       - inputqrs{4} = SV mapped detections
%  opt_sqi         Parameters calculated by previous function
%   fuse           If true good fidual candidates are fused using kernel
%                  density estimation [bool]
%
%% Parameters
WIN = 15;     % window used for segmenting signal previous to detection [sec]
THR = 0.5;     % threshold used by detector

%% Input test
if size(data,1) < size(data,2)
    data = data';
end

inputqrs(cellfun(@(x) isempty(x),inputqrs)) = [];
% accepted signals (pretty much everything that might make sense)
ecglike = {'eeg','eog','emg'};      % signals to be used
% bplike = {'pap' 'cvp'...};

%% Run

% figuring out what is what
idx = zeros(length(ecglike),length(header));
for i = 1:length(ecglike)
    idx(i,:) = cellfun(@any, regexpi(header,ecglike{i}));
end
ecgishleads = logical(sum(idx));

%=== testing for new available signals
if ~any(ecgishleads)
    qrsish = {};
else
    fprintf('\tEvaluating alternative signals... ');
    qrsish = cell(1,length(header));
    for k = find(ecgishleads)
        qrsish{k} = run_qrsdet_by_seg_fer(data(:,k),fs,WIN,THR);
%         figure
%         plot(data(:,k))
%         hold on
%         plot(qrsish{k},0.05,'or')
%         plot(data(:,1).*0.1,'g')
    end
end

%== ALL available QRSs
qrs_in = [inputqrs{:} qrsish{:}];
qrs_final = cell(opt_sqi.N_WIN,1);

%==== Regularity
M = numel(qrs_in);
% regularity!!
for w=1:opt_sqi.N_WIN
    curr_qrs = cell(1,M);
    ww = w*opt_sqi.REG_WIN; % in seconds
    SMI = 100*ones(1,numel(curr_qrs)+1);
    for m=1:M
        curr_qrs{m} = qrs_in{m}(qrs_in{m}>ww-opt_sqi.REG_WIN & qrs_in{m}<=ww);
        
        if numel(curr_qrs{m})>2
            SMI(m) = assess_regularity(curr_qrs{m}-(ww-opt_sqi.REG_WIN),0,1,0.96,opt_sqi.REG_WIN,0);
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


plot(data(:,1))
hold on
plot(qrs_final,1,'or')
close all

% 
% if fuse
%     qrs_final = kde_fusion(qrs_cand,fs,size(data,1));
% end



% det=kde_fusion(qrs)
%
% This functions uses kernel density estimation to fuse detections from
% different channels (or detectors).
%
% Input:
%  qrs              array with  samplestam for all detections (sorting
%                   is not relevant)
%  fs               signal sampling frequency
%  dlength          data length (peaks are bounded by size of measurement)
% Output:
%  det              fusioned detections
%
function det=kde_fusion(qrs,fs,dlength)
w_std = 0.05*fs;    % standard deviation of gaussian kernels [ms]
pt_kernel = round(fs/2);    % half window for kernel function [samples]

%% Calculating Kernel Density Estimation
% preparing annotations
peaks = hist(qrs,1:dlength);

% kde (adding gaussian kernels around detections)
kernel = exp(-(-pt_kernel:pt_kernel).^2./(2*w_std^2));

% calculating kernel density estimate
kde = conv(peaks,kernel,'same');

%% Decision
% Parameters
min_dist = round(0.4*fs);    % minimal distance between consecutive peaks [ms]
th = max(kde)/3;      % threshold for true positives (good measure is
% number_of_channels/3)

% Finding candidate peaks
wpoints = 1+min_dist:min_dist:dlength; % start points of windows (50% overlap)
wpoints(wpoints>dlength-2*min_dist) = []; % cutting edges
M = arrayfun(@(x) kde(x:x+2*min_dist-1)',wpoints(1:end), ...
    'UniformOutput',false);  % windowing signal (50% overlap)
M = cell2mat(M);
% adding first segment
head = [kde(1:min_dist) zeros(1,min_dist)]';
M = [head M];
% adding last segment
tail = kde((wpoints(end)+2*min_dist):dlength)';
tail = [tail; zeros(2*min_dist-length(tail),1)];
M = [M tail];
[~,idx] = max(M);   % finding maxima
idx = idx+[1 wpoints wpoints(end)+2*min_dist]; % absolute locations
i = 1;
while i<=length(idx)
    doubled = (abs(idx-idx(i))<min_dist);
    if sum(doubled) > 1
        [~,idxmax]=max(kde(idx(doubled)));
        idxmax = idx(idxmax + find(doubled,1,'first') - 1);
        idx(doubled) = [];
        idx = [idxmax idx];
        clear doubled idxmax
    end
    i = i+1;
end

% threshold check
idx(kde(idx)<th) = [];
det = sort(idx);
