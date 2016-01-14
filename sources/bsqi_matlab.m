function [ F1, IndMatch ] = bsqi_matlab( refqrs, testqrs,thres,fs)
%BSQI_MATLAB Calculate bSQI of two inputs
%   Detailed explanation goes here
%
% [1] Behar, Joachim, et al. "ECG signal quality during arrhythmia and its 
% application to false alarm reduction." Biomedical Engineering, IEEE Transactions on 60.6 (2013): 1660-1666.
%
% [2] Li, Qiao, Roger G. Mark, and Gari D. Clifford. "Robust heart rate estimation from multiple asynchronous noisy 
% sources using signal quality indices and a Kalman filter." Physiological measurement 29.1 (2008): 15.

% == managing inputs
if nargin<3; thres=0.05; end;
if nargin<4; fs=250; end;
thres = thres * fs;

    NB_REF  = length(refqrs);
    NB_TEST = length(testqrs);
% == core function
[IndMatch,Dist] = dsearchn(refqrs,testqrs);         % closest ref for each point in test qrs
IndMatchInWindow = IndMatch(Dist<thres);         % keep only the ones within a certain window
NB_MATCH_UNIQUE = length(unique(IndMatchInWindow)); % how many unique matching
TP = NB_MATCH_UNIQUE;                               % number of identified ref QRS
FN = NB_REF-TP;                                     % number of missed ref QRS
FP = NB_TEST-TP;                                    % how many extra detection?
Se  = TP/(TP+FN);
PPV = TP/(FP+TP);
F1 = 2*Se*PPV/(Se+PPV);                             % accuracy measure


end

