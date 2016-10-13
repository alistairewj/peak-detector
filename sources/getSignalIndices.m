function [ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header)
% [ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header)
%
% Given a 1xD cell array of strings, this function outputs the location of
% known signal types. The outputs returned are logical indices of the same
% size (1xD).
%
% It's essentially a set of regular expressions/string matches used to
% identify electrocardiogram signals from lead labels, for example, 'II'
% is usually Lead II of the ECG, and so idxECG is set to 1 for these
% elements, and so on.

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
ecg_leads = {'I','II','III','AVR','AVL',...
'AVF','V','V1','V2','V3','V4',...
'V5','V6','MCL','MCL1','MCL2','MCL3',...
'MCL4','MCL5','MCL6','aVR','aVL','aVF'};
idxECG = cellfun(@any, regexpi(header,'ecg','once'));
idxECG = idxECG | cellfun(@any, regexpi(header,'ekg','once'));

% exact match
idxECG = idxECG | ismember(header, ecg_leads);
% lowercase letters
idxECG = idxECG | ismember(header, lower(ecg_leads));
% uppercase letters
idxECG = idxECG | ismember(header, upper(ecg_leads));

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
