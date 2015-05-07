function [ann_sv, sq,rrDiff] = mapSVtoRR(ann_sv, ann_qrs, sq)
%=== find a minute long segment with a good matching
try
    N_START = 0; % seconds
    N_WINDOW = 60; % seconds
    %JOs
    if (nargin<3)
        sq = [];
    end
    
    if isempty(ann_qrs)
        ann_sv = ann_sv-0.2;
        idxRem = ann_sv<=0;
        ann_sv(idxRem) = [];
        if ~isempty(sq)
            sq(idxRem,:)=[];
        end
        rrDiff = 0.3;
        return;
    end
    foundMatch=false;
    while N_START < max(ann_qrs) && foundMatch == false
        N_END = N_START + N_WINDOW;
        idxKeep = ann_qrs>(N_START) & ann_qrs<(N_END);
        ann_qrs_map = ann_qrs(idxKeep);
        
        idxKeep = ann_sv>(N_START) & ann_sv<(N_END);
        ann_abp_map = ann_sv(idxKeep);
        
        if isempty(ann_abp_map) || isempty(ann_qrs_map)
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
        idxGood(idxGood>1) = 0; % treat these as bad matches
        if mean(idxGood)>0.9
            foundMatch = true;
            % this works, a lot of WABP match GQRS
        else
            N_START = N_START + N_WINDOW;
        end
    end
    
    %% map ABP to QRS
    if foundMatch==false
        rrDiff = 0.3;  % 200ms as a rough guess
    else
        ann_qrs_map = ann_qrs_map(idxMatch); % maps QRS to matching WABPs to the right of it
        [ann_qrs_map,idxUniq] = unique(ann_qrs_map);
        ann_abp_map = ann_abp_map(idxUniq); % remove double WABP detections if present
        rrDiff = mean(ann_abp_map-ann_qrs_map);
        
    end
    
    ann_sv = ann_sv - rrDiff;
    ann_sv = ann_sv(:);
    idxRem = ann_sv<=0;
    ann_sv(idxRem) = [];
    
    if nargin>2
        sq(idxRem,:) = [];
    end
catch
    ann_sv = {};
    sq   =  {};
    rrDiff = {};
end
