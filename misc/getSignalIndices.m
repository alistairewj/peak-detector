function [ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header)
% [ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header)
% looks for variations of signal names corresponding to one type of signal

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

