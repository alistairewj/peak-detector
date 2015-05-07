function [ sq, header, sqi, header_sqi ] = calcABPSQI(data, beat, fs)
%CALCABPSQI	Calculate ABP SQIs
%	[ sq, header ] = calcABPSQI(data, beat, samp_freq) calculates the
%   signal quality of the detected beats in the given waveform.
%
%	[ sq, header, sqi ] = calcABPSQI(data, beat, samp_freq) also outputs 
%	the SQIs calculated for the given beats.
%
%	[ sq, header, sqi, header_sqi ] = calcABPSQI(data, beat, samp_freq) 
%	also outputs a header cell array of strings for the sqi output
%
%	Inputs:
%       data		- Nx1 vector of ABP signal
%       beat		- Qx1 vector of samples of detected beats
%       fs 			- Sampling frequency
%
%	Outputs:
%       sq			- QxD vector of good or bad quality flags
%       header		- 1xD cell array of strings describing signal quality
%       sqi			- Qx(D+2) matrix of signal quality features (+time)
%       header		- 1xD+2 cell array of strings describing features
%
%	See also JSQI ABPFEATURE

%	Copyright 2012 Alistair Johnson

% References:
%   This code is based off of the SQIs developed in the following
%   publications.
%   Sun JX. Cardiac Output Estimation using Arterial Blood Pressure Waveforms. [MEng thesis]. Cambridge, MA: Massachusetts Institute of Technology, Department of Electrical Engineering and Computer Science, September 2006.
%   Sun JX, Reisner AT, Saeed M, Mark RG. Estimating Cardiac Output from Arterial Blood Pressure Waveforms: a Critical Evaluation using the MIMIC II Database. Computers in Cardiology 32:295-298 (2005).

header_sqi = {'Time of systole   [time]',...
'Systolic BP       [mmHg]',...
'Time of diastole  [time]',...
'Diastolic BP      [mmHg]',...
'Pulse pressure    [mmHg]',...
'Mean pressure     [mmHg]',...
'Beat Period       [time]',...
'mean_dyneg',...
'End of systole time  0.3*sqrt(RR)  method',...
'Area under systole   0.3*sqrt(RR)  method',...
'End of systole time  1st min-slope method',...
'Area under systole   1st min-slope method',...
'Pulse              [time]'};


    header={'Combined','P','MAP','HR','PP','Psys','Pdias','period','Ponset','noisy'};
%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 20-Jul-2012 12:36:37
%	Contact: alistairewj@gmail.com

sqi = abpfeature(data,beat);
if isempty(sqi) % not enough data
    sqi = zeros(size(beat,1),numel(header_sqi));
    sq = zeros(size(beat,1),numel(header));
    return;
end

sqi = [sqi;zeros(1,size(sqi,2))]; % last beat 0 signal quality
%           Col 1:  Time of systole   [samples]
%               2:  Systolic BP       [mmHg]
%               3:  Time of diastole  [samples]
%               4:  Diastolic BP      [mmHg]
%               5:  Pulse pressure    [mmHg]
%               6:  Mean pressure     [mmHg]
%               7:  Beat Period       [samples]
%               8:  mean_dyneg
%               9:  End of systole time  0.3*sqrt(RR)  method
%              10:  Area under systole   0.3*sqrt(RR)  method
%              11:  End of systole time  1st min-slope method
%              12:  Area under systole   1st min-slope method
%              13:  Pulse             [samples]
[sq, r] = jSQI(sqi, beat, data, fs);

if size(sqi,1)~=size(beat,1)
    fprintf('sqi: %d\tbeat: %d\n',size(sqi,1),size(beat,1));
end
if size(sq,1)~=(size(beat,1))
    fprintf('sqi: %d\tbeat: %d\n',size(sq,1),size(beat,1));
end
% update SQI values in samples to time, using fs which should be 125 Hz due
% to resampling
sqi(:,[1,3,7,9,11]) = sqi(:,[1,3,7,9,11]) / fs;
sqi(:,13) = sqi(:,13) * fs;
sqi = [sqi;zeros(1,size(sqi,2))];

end