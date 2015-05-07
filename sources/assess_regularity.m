function [SMI,hr] = assess_regularity(QRS,secDer,fs,CI,segL,debug)
% compute smoothness of hrv or of hr time series given a confidence interval (CI).
% The underlying assumption for using this function is that the more smooth
% the QRS time series is the most likely it is to be a meaningful QRS time series. 
%
% inputs
%   QRS:    QRS fiducials (required, in data points number)
%   secDer: use second derivative? (i.e. HRV instead of HR to 
%           compute SMI, default: 1)
%   fs:     sampling frequency (default: 1kHz)
%   CI:     confidence interval (default: 0.96)
%   segL:   length of the ecg segment in seconds (default: 60sec)
%
% outputs
%   SMI:    regularity measure which corresponds to the standard deviation 
%           of hrv OR hr (depending on secDer option)
%   hr:     heart rate
%
%
% Safe Foetus Monitoring Toolbox, version 1.0, Sept 2013
% Released under the GNU General Public License
%
% Copyright (C) 2013  Joachim Behar
% Oxford university, Intelligent Patient Monitoring Group - Oxford 2013
% joachim.behar@eng.ox.ac.uk
%
% Last updated : 02-08-2013
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.

% == manage inputs
if nargin<1; error('assess_regularity: wrong number of input arguments \n'); end;
if nargin<2; secDer = 1; end;
if nargin<3; fs = 1000; end;
if nargin<4; CI = 1; end;
if nargin<5; segL = 60; end;
if nargin<6; debug = 0; end;

% == core function
%try
    % == compute variability
    hr = 60./(diff(QRS)/fs);

    if secDer
        % rather than looking at the distribution of hr this option is looking at
        % the distribution of hrv. This makes more sense because this way
        % we are looking into high variation in deltaHR from a measure to the
        % following one rather than the variability of absolute value of HR (which 
        % might be high if the foetus HR is changing significantly)

        % use the second derivative of QRS location i.e. ~diff(HR) i.e. hrv
        % because we do not know how many QRS have been detected we are going
        % to resample at a defined number of datapoints xi
        xi = fs*(0.1:0.45:segL-1);
        % hr is resampled using a spline function to smooth it. So the output 
        % is a smoothed version of the hr (yi) at a defined number of
        % datapoints (xi)
        yi = interp1(QRS(1:end-1),hr,xi,'spline');
        % now taking the derivative of smoothed hear rate. This will give the
        % hrv
        hrv = sort(diff(yi));
        hrv_N = length(hrv);
        % plot(hr); hold on, plot(yi,'r');
        % we tolerate some mistakes using a confidence interval
        if CI~=1
            CI_sup = ceil(hrv_N*(CI+(1-CI)/2));
            CI_inf = ceil(hrv_N*((1-CI)/2));
            hrv_CI = hrv(CI_inf:CI_sup);
        else
            hrv_CI = hrv;
        end
        % output the std of the hrv 
        % SMI = std(hrv_CI); % OLD version

        % new version (25-08-2013)
        SMI = length(find(hrv_CI>30 | hrv_CI<-30));
        % this looks at the absolute number of outliers in hrv_CI with >
        % 30bpm drop or increase from a point to the next.

    else
        % use first derivative of RR time series i.e. hr
        hr_sorted = sort(hr);
        hr_N  = length(hr_sorted);
	if hr_N < 3
		SMI = 100;
	else
        	if CI~=1
            CI_sup = ceil(hr_N*(CI+(1-CI)/2));
            CI_inf = ceil(hr_N*((1-CI)/2));
            hr_CI  = hr_sorted(CI_inf:CI_sup);
        	else
            hr_CI  = hr_sorted;
        	end
        SMI = std(hr_CI);
	end
end

%catch ME
%    for enb=1:length(ME.stack); disp(ME.stack(enb)); end;
%    SMI = 100; hr = [];
%end

% == plots
if debug
    hist(hrv_CI,40); xlabel('hr or hrv histogram');
    title(['assess regularity in term of hr or hrv, REGULARITY:' SMI]);
    set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold');
end

end

