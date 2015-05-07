function [ann_abp, sq, rrDiff, foundMatch] = mapABPtoRR(ann_abp, ann_qrs,sq)
%=== find a minute long segment with a good matching
N_START = 0; % seconds
N_WINDOW = 60; % seconds

if nargin<3
    sq=[];
end
%JOs 25/02/2014 move initialisation fo foudnMatch in case ann_qrs is empty
foundMatch=false;
if isempty(ann_qrs)
    ann_abp = ann_abp-0.2;
    idxRem = ann_abp<=0;
    ann_abp(idxRem) = [];
    if nargin>2
        sq(idxRem,:)=[];
    end
    rrDiff = 0.2;
    return;
end

while N_START < max(ann_qrs) && foundMatch == false
    %%
    N_END = N_START + N_WINDOW;
    idxKeep = ann_qrs>(N_START) & ann_qrs<(N_END);
    ann_qrs_map = ann_qrs(idxKeep);
    
    idxKeep = ann_abp>(N_START) & ann_abp<(N_END);
    ann_abp_map = ann_abp(idxKeep);
    
    if numel(ann_abp_map)<10 || numel(ann_qrs_map)<10
        N_START = N_START + N_WINDOW;
        continue;
    end
    
    % no QRS annotations within ABP annotations
    if all(ann_abp_map < ann_qrs_map(1) | ann_abp_map > ann_qrs_map(end))
        N_START = N_START + N_WINDOW;
        continue;
    end
    
    
    % make sure ann_gqrs circumscribes wabp
    while ann_abp_map(1)<ann_qrs_map(1)
        ann_abp_map(1) = [];
    end
    
    while ann_abp_map(end)>ann_qrs_map(end)
        ann_abp_map(end) = [];
    end
    
    [idxGood,idxMatch] = histc(ann_abp_map,ann_qrs_map);
    idxGood(idxGood>1) = 0;
    
    if mean(idxGood)>0.9 % at least 10 good matches
        foundMatch = true;
        % this works, a lot of WABP match GQRS
    else
        N_START = N_START + N_WINDOW;
    end
end

%% map ABP to QRS
if foundMatch==false
    rrDiff = 0.2; % 200ms as a rough guess
else
    ann_qrs_map = ann_qrs_map(idxMatch); % maps QRS to matching WABPs to the right of it
    [ann_qrs_map,idxUniq] = unique(ann_qrs_map);
    ann_abp_map = ann_abp_map(idxUniq); % remove double WABP detections if present
    rrDiff = median(ann_abp_map-ann_qrs_map);
end

ann_abp = ann_abp - rrDiff;
ann_abp = ann_abp(:);
idxRem = ann_abp<=0;
ann_abp(idxRem) = [];

if nargin>2 && ~isempty(sq)
   sq(idxRem,:) = [];
end

end
