function QRS = run_qrsdet_by_seg_ali(ecg,fs,opt)
% this function is used to run the QRS detector for each window window (non overlapping)
% this is needed because in the case of big artefacts at the begining of
% the record (i.e. where the P&T threshold is computed) it can make the detection fail.
%
% inputs
%   ecg:        ecg signal
%   fs:         sampling frequency
%   window:     size of the window onto which to perform QRS detection (in seconds)
%   thres:      threshold to be used
% 
% output
%   QRS:        QRS location in nb samples (ms)
%
%
% FECG-ESN toolbox, version 1.0
% Released under the GNU General Public License
%
% Copyright (C) 2014  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 28-01-2014
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


% == managing inputs
if nargin<1; error('run_qrsdet_by_seg: wrong number of input arguments \n'); end;
if nargin<2; fs=1000; end;
if nargin<3 || ~isstruct(opt)
    opt = setOptions;
else
    [ opt ] = setOptions(opt);
end

% == general
segsizeSamp = opt.JQRS_WINDOW*fs; % convert window into nb of samples
NbSeg = floor(length(ecg)/segsizeSamp); % nb of segments
QRS = cell(NbSeg,1);
signForce = 0; % if we want to force the sign of the peak we are looking for



if NbSeg == 0
    QRS = 1;
    return;
end
    
% == core function
try
    
    %=== First subsegment
    % first subsegment - look forward 1s
    dTplus  = fs;
    dTminus = 0;
    start = 1;
    stop  = segsizeSamp;
    
    %=== if no more data, don't look ahead
    if NbSeg==1
        dTplus  = 0;
        stop  = length(ecg);
    end
    
    % sign of peaks is determined by the sign on the
    % first window and then is forced for the following windows.
    [QRStemp,signForce] = qrs_detect2(ecg(start-dTminus:stop+dTplus),opt.JQRS_REFRAC,opt.JQRS_THRESH,fs,[],signForce,0,opt.JQRS_INTWIN_SZ);
    QRS{1} = QRStemp(:);
    
    start = start+segsizeSamp;
    stop = stop+segsizeSamp;
    
    % for each segment perform QRS detection
    for ch=2:NbSeg-1

        % take +/-1sec around selected subsegment exept for the borders. This
        % is in case there is a QRS in between segments -> allows to locate
        % them well.
        dTplus  = fs;
        dTminus = fs;
        
        [QRStemp,signForce] = qrs_detect2(ecg(start-dTminus:stop+dTplus),opt.JQRS_REFRAC,opt.JQRS_THRESH,fs,[],signForce);

        NewQRS = (start-1)-dTminus+QRStemp;
        NewQRS(NewQRS>stop) = [];
        NewQRS(NewQRS<start) = [];

        if ~isempty(NewQRS) && ~isempty(QRS{ch-1})
            % this is needed to avoid multiple detection at the transition point
            NewQRS(NewQRS<QRS{ch-1}(end)) = [];
            if ~isempty(NewQRS) && (NewQRS(1)-QRS{ch-1}(end))<opt.JQRS_REFRAC*fs
                % between two windows
                NewQRS(1) = [];
            end
            
        end
        QRS{ch} = NewQRS(:);

        start = start+segsizeSamp;
        stop = stop+segsizeSamp;
    end
    
    %JOs 24/02/15 check there is more than one segment otherwise
    %PROBLEM...
    if NbSeg>1
        % last subsegment
        ch = NbSeg;
        stop  = length(ecg);
        dTplus  = 0;
        dTminus = fs;
        [QRStemp,signForce] = qrs_detect2(ecg(start-dTminus:stop+dTplus),opt.JQRS_REFRAC,opt.JQRS_THRESH,fs,[],signForce);
        
        NewQRS = (start-1)-dTminus+QRStemp;
        NewQRS(NewQRS>stop) = [];
        NewQRS(NewQRS<start) = [];
        
        if ~isempty(NewQRS) && ~isempty(QRS{ch-1})
            % this is needed to avoid multiple detection at the transition point
            NewQRS(NewQRS<QRS{ch-1}(end)) = [];
            if ~isempty(NewQRS) && (NewQRS(1)-QRS{ch-1}(end))<opt.JQRS_REFRAC*fs
                % between two windows
                NewQRS(1) = [];
            end
            
        end
        QRS{ch} = NewQRS(:);
    end
    
    %=== convert to double
    QRS = vertcat(QRS{:});

catch ME
    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
%    QRS = [1000 2000];
    rethrow(ME);
end

end




function [ opt ] = setOptions(opt)
opt_default=struct();
opt_default.JQRS_THRESH = 0.3;
opt_default.JQRS_WINDOW = 15;
opt_default.JQRS_REFRAC = 0.250;
opt_default.JQRS_INTWIN_SZ = 7;
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
