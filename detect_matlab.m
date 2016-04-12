function [ qrs, sqi ] = detect_matlab(data, header, fs, alg, opt)
%[ qrs, sqi ] = detect_matlab(data, header, fs, alg, opt) detects QRS 
%complexes in the provided data. The QRS detections and signal quality
%estimates are ouput.

% Dependencies:
%
%       1) This function requires the WFDB Toolbox for MATLAB and Octave.
%          For information on how to install the toolbox, please see:
%              http://www.physionet.org/physiotools/matlab/wfdb-app-matlab/
%
%       2) The Toolbox is supported only on 64-bit MATLAB 2013a and 2013b
%          and on GNU Octave 3.6.4 and later, on Linux, Mac OS X, and
%          Windows.

%% CHECK INPUTS ARE VALID

if nargin<3
    error('Function requires at least 3 inputs.');
end

% check data types
if isnumeric(data)==0
    error('detect_matlab:invalidData',...
        'Data must be an NxD matrix of N samples and D signals.');
end
if ischar(header) && size(data,2)==1
    header = {header}; % convert to 1x1 cell array of strings
elseif ~iscell(header)
    error('detect_matlab:invalidHeader',...
        'Header must be an 1xD cell array of strings.');
end

% check that data/header are consistent in size
if size(data,2) ~= size(header,2)
    error('detect_matlab:dataMismatch',...
        'Each signal in DATA must have a corresponding element in HEADER, and vice versa.');
end

if isnumeric(fs)==0
    error('detect_matlab:invalidSamplingFrequency',...
        'Third input must be a scalar or vector of integers representing the sampling frequency.');
end
if numel(fs) == 1
    fs = repmat(fs,1,size(data,2));
end

if nargin<4
    alg = 'sqi';
elseif ~ischar(alg)
    error('detect:incorrectAlgorithm',...
        'Switching algorithm should be a string: either ''sqi'' or ''regularity''.')
elseif ~strcmpi(alg,'sqi') && ~strcmpi(alg,'regularity')
    error('detect:incorrectAlgorithm',...
        'Switching algorithm should be ''sqi'' or ''regularity''')
end

% if no options input, create empty structure
% option settings are checked in the detect subfunctions
if nargin<5
    opt = struct();
end

%% run either SQI algorithm or regularity algorithm
if strcmpi(alg,'regularity')
    % run heart rate regularity based switching
    [ qrs, sqi, qrs_comp, qrs_header ] = detect_regularity(data, header, fs, opt);
else
    % run SQI based switching (it's the only other option)
    [ qrs, sqi, qrs_comp, qrs_header ] = detect_sqi(data, header, fs, opt);
end

end
