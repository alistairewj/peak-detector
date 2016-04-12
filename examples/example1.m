% Example showing how to use the detection function with WFDB records
recordName = '100_snip';

% Note, by default, the WFDB toolbox looks in the local directory for data
% So we create a variable to point it to the data
% If you do not run this from the base folder (usually peak-detector), you
% will need to update this path.
data_path = 'examples/';

% define input options for the peak detector
% all of the options listed here are the default values, and are optionally omitted
opt = struct(...
    'SIZE_WIND',10,... % define the window for the bSQI check on the ECG
    'LG_MED',3,... % take the median SQI using X nearby values,  so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows
    'REG_WIN',1,... % how frequently to check the SQI for switching - i.e., if REG_WIN = 1, then we check the signals every second to switch
    'THR',0.150,... % the width, in seconds, used when comparing peaks in the F1 based ECG SQI
    'SQI_THR',0.8,... % the SQI threshold - we switch signals if SQI < this value
    'USE_PACING',1,... % flag turning on/off the pacing detection/correction
    'ABPMethod','wabp',... % ABP peak detection method (wabp, delineator)
    'SIMPLEMODE',0,... % simple mode only uses the first ABP and ECG signal, and ignores all others
    'DELAYALG', 'map',... % algorithm used to determine the delay between the ABP and the ECG signal
    ... % jqrs parameters - the custom peak detector implemented herein
    'JQRS_THRESH', 0.3,... % energy threshold for defining peaks
    'JQRS_REFRAC', 0.25,... % refractory period in seconds
    'JQRS_INTWIN_SZ', 7,...
    'JQRS_WINDOW', 15);

[ qrs, sqi ] = detect([data_path recordName], 'regularity', opt);
