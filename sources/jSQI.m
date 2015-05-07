function [BeatQ, r] = jSQI(features, onset, abp, fs)
% JSQI  ABP waveform signal quality index.
%   [BEATQ, R] = JSQI(FEATURES, ONSET, ABP,fs) returns a binary signal quality
%   assessement of each beat in ABP.  This algorithm relies on detecting
%   abnormalities of numeric values in FEATURES and ONSET.
%
%   In:   FEATURES <mx12> numeric matrix --- features extracted from ABP using abpfeature.m
%         ONSET    <nx1>  numeric vector --- onset times of ABP using wabp.m
%         ABP      <px1>  numeric vector --- arterial blood pressure waveform
%         FS       <1x1>  numeric double --- sampling frequency
%
%   Out:  BEATQ    <nx10> logical matrix --- SQI of each beat: 0=good, 1=bad
%           Col 1:  logical OR of cols 2 thru 10
%               2:  P   not physiologic  (<20 or >300 mmHg)
%               3:  MAP not physiologic  (<30 or >200 mmHg)
%               4:  HR  not physiologic  (<20 or >200 bpm)
%               5:  PP  not physiologic  (<30 mmHg)
%               6:  abnormal Psys        (beat-to-beat change > 20 mmHg)
%               7:  abnormal Pdias       (beat-to-beat change > 20 mmHg)
%               8:  abnormal period      (beat-to-beat change > 1/2 sec)
%               9:  abnormal P(onset)    (beat-to-beat change > 20 mmHg) 
%              10:  noisy beat           (mean of negative dP < -3)
%         R        <1x1>  fraction of good beats in ABP
% 
%   Usage:
%   - FEATURES must be obtained using abpfeature.m
%   - ONSET    must be obtained using wabp.m
%
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.
%   - v2.0 - 1/18/06 - thresholds updated to reduce false positives
%   - v3.0 - 2/10/06 - added "..101..." detection - see lines 92-96
%
%   Updated by Alistair Johnson (alistairewj@gmail.com) on April 1st, 2014
%   - Added sampling frequency as input. This effects:
%       The HR calculation
%       The time based calculations (dPOnset, dPeriod)
%   - Increased heart rate threshold to 240 bpm.
%   - Modified the overall signal quality calculation to use a subset of the features

if nargin<4
    fs = 125;
end

if length(onset)<30
    BeatQ = [];
    r = [];
    return
end

% thresholds
rangeP      = [20 300]; % mmHg
rangeMAP    = [30 200]; % mmHg
rangeHR     = [20 240]; % bpm
rangePP     = [30 inf]; % mmHg

dPsys       = 20;           % Maximum beat to beat change in systole
dPdias      = 20;           % Maximum beat to beat change in diastole
dPeriod     = 0.5*fs;       % 1/2 second
dPOnset     = 20;           % Maximum change in onset of pulse
%dPOnset - add in fs?
noise       = -3;


% get ABP features
Psys        = features(:,2);
Pdias       = features(:,4);
PP          = features(:,5);
MAP         = features(:,6);
BeatPeriod  = features(:,7);
mean_dyneg  = features(:,8);
HR          = features(:,13)*fs;

% absolute thresholding (flag unphysiologic beats)
badP        = find(Pdias < rangeP(1)   |   Psys > rangeP(2));
badMAP      = find(MAP   < rangeMAP(1) |    MAP > rangeMAP(2));
badHR       = find(HR    < rangeHR(1)  |     HR > rangeHR(2));
badPP       = find(PP    < rangePP(1));

% first difference thresholding (flag beat-to-beat variations)
jerkPsys    = 1 + find(abs(diff(Psys))       > dPsys);
jerkPdias   =     find(abs(diff(Pdias))      > dPdias);
jerkPeriod  = 1 + find(abs(diff(BeatPeriod)) > dPeriod);
jerkPOnset  =     find(abs(diff(abp(onset))) > dPOnset);

% noise detector
noisy       = find(mean_dyneg < noise);

% SQI final
bq = zeros(length(onset),10);
bq(badP,       2) = 1;
bq(badMAP,     3) = 1;
bq(badHR,      4) = 1;
bq(badPP,      5) = 1;
bq(jerkPsys,   6) = 1;
bq(jerkPdias,  7) = 1;
bq(jerkPeriod, 8) = 1;
bq(jerkPOnset, 9) = 1;
bq(noisy,     10) = 1;

bq(:,1) = any(bq(:,2:5),2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % make all "...101..." into "...111..." (aka, smooth out transients)
% y = bq(:,1);
% y(find(diff(y,2)==2)+1)=1;
% bq(:,1)=y;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BeatQ   = bq==0; % 1=good signal quality
% fraction of good beats overall
r = length(find(bq(:,1)==1))/length(onset);

end
