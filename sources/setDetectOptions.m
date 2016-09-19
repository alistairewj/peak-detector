function [ opt ] = setDetectOptions(varargin)
% Creates a structure of options used in the SQI/peak detection algorithm
% Input should be parameter value pairs
% If the string does not match an option, that parameter value pair is
% ignored
%
% Example:
%   opt = setDetectOptions('LG_MED',3,'REG_WIN',1);

%=== parameters that define the window for the bSQI check on the ECG
opt_default.SIZE_WIND = 10;

% take the median SQI using X nearby values
% this is used only for the ECG SQI
opt_default.LG_MED = 3;
% so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows

% how frequently to check the SQI for switching
% i.e., if REG_WIN = 1, then we check the signals every second to switch
opt_default.REG_WIN = 1;

% the width, in seconds, used when comparing peaks in the F1 based ECG SQI
opt_default.THR = 0.150;

% the SQI threshold
%   if the lead SQI is higher than this we use this signal
%   if the lead SQI is lower than this, we use the next signal
%   this an ordered comparison, so we usually default to the ECG signal
%   (the first column/signal present in the data)
opt_default.SQI_THR = 0.8;
opt_default.USE_PACING = 1; % flag turning on/off the pacing detection/correction

% ABP peak detection method
%   options:
%       wabp
%       delineator
opt_default.ABPMethod = 'wabp';

opt_default.SIMPLEMODE = 0; % only use the first of ABP/ECG

%=== jqrs parameters
opt_default.JQRS_THRESH = 0.3;
opt_default.JQRS_REFRAC = 0.25;
opt_default.JQRS_INTWIN_SZ = 7;
opt_default.JQRS_WINDOW = 15;

opt_default.DELAYALG = 'map';

% Flag which turns off SV/PPG detectors
% these detectors are currently not used in the 'sqi' approach
opt_default.ENABLE_OTHER_DETECTORS = 0;

% leave temporary files in the working directory
opt_default.SAVE_STUFF = 0;

if nargin==0
    opt = opt_default;
    return;
elseif nargin==1
    opt = varargin{1};
else
    error('Incorrect number of inputs.');
end

% if input options are given, we update the default opt structure with them
% we also do a lot of input argument checking
if nargin>0 && isstruct(opt)
    fn = fieldnames(opt);
    fn_default = fieldnames(opt_default);
    for f=1:numel(fn)
        if ismember(fn{f},fn_default)
            val = opt.(fn{f});
            
            switch fn{f}
                case {'USE_PACING','SAVE_STUFF','SIMPLEMODE','ENABLE_OTHER_DETECTORS'} 
                    % binary
                    if ~isnumeric(val) || (val ~= 0 && val ~=1)
                        error('setDetectOptions:badValue',...
                            '%s: field can only take values 0 or 1.',...
                            fn{f});
                    end
                case {'THR','JQRS_THRESH','JQRS_REFRAC'} 
                    % numeric
                    if ~isnumeric(val)
                        error('setDetectOptions:badValue',...
                            '%s: field should be an integer.',...
                            fn{f});
                    end
                        
                case {'SQI_THR',} 
                    % numeric between 0-1
                    if ~isnumeric(val) || val < 0 || val > 1
                        error('setDetectOptions:badValue',...
                            '%s: field should be an integer.',...
                            fn{f});
                    end
                    
                
                case {'LG_MED','REG_WIN','JQRS_INTWIN_SZ','JQRS_WINDOW'}
                    % ensure integer
                    if ~isnumeric(val) || round(val) ~= val
                        error('setDetectOptions:badValue',...
                            '%s: field should be an integer.',...
                            fn{f});
                    end
                    
                case {'SIZE_WIND'}
                    % ensure integer *AND* divisible by 2
                    if ~isnumeric(val) || round(val) ~= val ||  mod(val,2) ~= 0
                        error('setDetectOptions:badValue',...
                            '%s: field should be an integer divisible by two.',...
                            fn{f});
                    end
                case {'ABPMethod','DELAYALG'}
                    % special case strings
                    if ~ischar(val)
                        error('setDetectOptions:badValue',...
                            '%s: field should be a string.',...
                            fn{f});
                    end
                    
                    switch fn{f} 
                        case 'ABPMethod'
                            if any(ismember(val,{'wabp','delineator'}))==0
                                error('setDetectOptions:badValue',...
                                    '%s: field must be either ''wabp'' or ''delineator''.',...
                                    fn{f});
                            end
                        case 'DELAYALG'
                            if any(ismember(val,{'map','crosscorr'}))==0
                                error('setDetectOptions:badValue',...
                                    '%s: field must be either ''map'' or ''crosscorr''.',...
                                    fn{f});
                            end
                    end
                                
                    
            end
            
            opt_default.(fn{f}) = val;
        else
            fprintf('Ignoring unrecognized option %s\n',fn{f});
        end
    end
end

% some options are derived from others, and used for convenience
opt_default.HALF_WIND = opt_default.SIZE_WIND/2;
opt = opt_default;

end