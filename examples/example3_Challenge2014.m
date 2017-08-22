% This example runs the peak detector on data from the PhysioNet/CinC 2014
% challenge. The aim of this challenge was multimodal peak detection.
% You will need to have downloaded and extracted the challenge data.
% You can find the challenge description here:
%    https://physionet.org/challenge/2014/#challenge-data
% You can download the training data directly using this link:
%   https://physionet.org/physiobank/database/challenge/2014/set-p.zip
% There is another 100 records from the challenge which you can also use:
%   https://physionet.org/challenge/2014/training.zip

% Requirements
%   You will need the WFDB for MATLAB toolbox.
%   Specifically, you need to add the mcode folder to your path.

% First, specify the path with the data
% The path must have a RECORDS file, which lists all records in the folder
challenge_data_path = '/data/challenge2014/set-p';



if exist(challenge_data_path, 'dir')~=7
    error('File path does not exist');
end

records_file = [challenge_data_path,filesep,'RECORDS'];
if exist(records_file, 'file')~=2
    error('Folder must contain RECORDS file.');
end

% HACK: WFDB_PATH is the default search path for WFDB.
% It can be manually changed in wfdbloadlib.m
% We use change directory to avoid dealing with this setting.
currentFolder = pwd;

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


% loop through each record
fp = fopen(records_file);
records = textscan(fp,'%s');
records = records{1};
fclose(fp);

cd(challenge_data_path);
figure(1);
for r = 1:numel(records);
    try
        % Example showing how to use the detection function with WFDB records
        recordName = records{r};
        fprintf('=====\n%s\n=====\n', recordName);


        % load the data into matlab
        data = rdsamp(recordName);
        [siginfo,fs,sigClass]=wfdbdesc(recordName);
        header = arrayfun(@(x) x.Description, siginfo, 'UniformOutput', false);

        % pass the data to the detect_sqi function 
        [ qrs, sqi, qrs_comp, qrs_header ] = detect_sqi(data, header, fs, opt);
        
        % NOTE: we could have called detect() with the record name as
        % follows:
        % [ qrs, sqi ] = detect(recordName, 'sqi', opt);
        % however, since we want to plot the data, we load above
        
        % Plot the output
        set(0,'CurrentFigure',1); clf;
        plotann( data, header, fs(1), qrs_comp, sqi );
        drawnow;
    catch ME
        cd(currentFolder);
        rethrow(ME);
    end
end
cd(currentFolder);



