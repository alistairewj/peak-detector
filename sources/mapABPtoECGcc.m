function [delay, ann_abp] = mapABPtoECGcc(ecg, abp, fs, ann_abp)
% Maps ABP to ECG using cross correlation
% Optionally updates ABP annotations with the delay

if nargin<4
    ann_abp = [];
end

WINDOW_SIZE = 20; % seconds
N = numel(ecg);
T = N ./ fs; % maximum time in seconds
N_LAG = 3*fs; % set number of lags equal to fs, so maximum lag is 1 second
N_WINDOW = ceil(T ./ WINDOW_SIZE);
DEBUG=0;
%=== let's preprocess the signal a little bit
[ ecg ] = preprocessTheSignalALittleBit(ecg, fs);

if isempty(ecg)
    fprintf('ECG is mostly a flat line. Providing a default delay of 200ms.\n');
    delay = ceil(ones(N_WINDOW,1)*fs*0.2);
    return;
else
    %=== First, calculate a distribution of delays
    delay = zeros(N_WINDOW,1);
end
w = 1;
while w<N_WINDOW
    idx = (w-1)*WINDOW_SIZE*fs+1 : w*WINDOW_SIZE*fs;
    [ amp, lag ] = mycrosscorr(ecg(idx),abp(idx),N_LAG);
    %     figure(2); clf;stem(lag,c);
    
    [~,idxLag] = max(amp);
    delay(w) = lag(idxLag);
    if idxLag==numel(amp)
        fprintf('Lag %d is at the end of the cross-corr. Implies >1 sec delay.\n',lag(idxLag));
    end
    %% debugplot
    if DEBUG==1
        figure(1); clf;
        col = [0.904705882352941,0.191764705882353,0.198823529411765;0.294117647058824,0.544705882352941,0.749411764705882;0.371764705882353,0.717647058823529,0.361176470588235;1,0.548235294117647,0.100000000000000;0.955000000000000,0.894647058823529,0.472176470588235;0.685882352941177,0.403529411764706,0.241176470588235;0.971764705882353,0.555294117647059,0.774117647058824];
        subplot(2,1,1);stem(lag,amp);subplot(2,1,2);
        hold all;
        plot((idx-idx(1))/fs,ecg(idx)./prctile(ecg,99),'-','linewidth',2,'color',col(1,:));
        plot((idx-idx(1))/fs,abp(idx)./prctile(abp,99),'-','linewidth',2,'color',col(2,:));
        plot((idx+delay(w)-idx(1))/fs,abp(idx)./prctile(abp,99),'-','linewidth',2,'color',col(3,:));
    end
    w=w+1;
end

%=== for the last window, we only use it if we have at least 5 sec data
if (N-((w-1)*WINDOW_SIZE*fs)) > (fs*5)
    idx = (w-1)*WINDOW_SIZE*fs+1 : N;
    [ amp, lag ] = mycrosscorr(ecg(idx),abp(idx),N_LAG);
    [~,idxLag] = max(amp);
    if idxLag~=numel(amp)
        delay(w) = lag(idxLag);
    else
        fprintf('Lag %d is at the end of the cross-corr. Implies >1 sec delay.\n',lag(idxLag));
    end
else
    delay = delay(1:end-1);
end

%=== delete the bad ones which are too long
delay(delay==lag(end)) = [];

% %=== distribution of the lags
% figure(3); clf; hist(delay,0:1:fs);

end

function [ ecg ] = preprocessTheSignalALittleBit(ecg, fs)

ecg = ecg(:); % ensure column vector
MIN_AMP = 0.1; % if the median of the filtered ECG is inferior to MINAMP then it is likely to be a flatline
% note the importance of the units here for the ECG (mV)
N = numel(ecg); % number of input samples
MED_SMOOTH_NB_COEFF = round(fs/100);
INT_NB_COEFF = round(7*fs/256); % length is 7 for fs=256Hz

% == Bandpass filtering for ECG signal
% this sombrero hat has shown to give slightly better results than a
% standard band-pass filter. Plot the frequency response to convince
% yourself of what it does
b1 = [-7.757327341237223e-05  -2.357742589814283e-04 -6.689305101192819e-04 -0.001770119249103 ...
    -0.004364327211358 -0.010013251577232 -0.021344241245400 -0.042182820580118 -0.077080889653194...
    -0.129740392318591 -0.200064921294891 -0.280328573340852 -0.352139052257134 -0.386867664739069 ...
    -0.351974030208595 -0.223363323458050 0 0.286427448595213 0.574058766243311 ...
    0.788100265785590 0.867325070584078 0.788100265785590 0.574058766243311 0.286427448595213 0 ...
    -0.223363323458050 -0.351974030208595 -0.386867664739069 -0.352139052257134...
    -0.280328573340852 -0.200064921294891 -0.129740392318591 -0.077080889653194 -0.042182820580118 ...
    -0.021344241245400 -0.010013251577232 -0.004364327211358 -0.001770119249103 -6.689305101192819e-04...
    -2.357742589814283e-04 -7.757327341237223e-05];

b1 = resample(b1,fs,250);
ecg = filtfilt(b1,1,ecg);

%     N_GOOD = sum(abs(ecg)>MIN_AMP);
%     if (N_GOOD/N)>0.20 || ...
%             N_GOOD > (fs*2*60)
% if 20% of the samples have an absolute amplitude which is higher
% than MIN_AMP then we are good to go.
%   OR
% if we have at least 2 minutes of data


if  abs( prctile(ecg,99) - prctile(ecg,1) ) >= MIN_AMP
    % we have some reasonable consistent variation in the signal
    
    % == P&T operations
    ecg = [ecg(2)-ecg(1);diff(ecg)];  % (4) differentiate (one datum shorter)
    ecg = ecg.*ecg; % (5) square ecg
    ecg = filtfilt(ones(1,INT_NB_COEFF),1,ecg); % (6) integrate
    ecg = medfilt1(ecg,MED_SMOOTH_NB_COEFF);  % (7) smooth
    
else
    ecg = []; % indicates the ECG is too poor quality to use
end
end

function [ amp, lag ] = mycrosscorr(ecg,abp,N_LAG)
%=== get FFT
N = numel(ecg);
NFFT = 2^(nextpow2(N));
f = fft([ecg-mean(ecg),abp-mean(abp)],NFFT);

%=== auto correlation = PSD(1)
autocorr_ecg = ifft(f(:,1).*conj(f(:,1)));
autocorr_ecg = sqrt(autocorr_ecg(1));
autocorr_abp = ifft(f(:,2).*conj(f(:,2)));
autocorr_abp = sqrt(autocorr_abp(1));

%=== multiply FFT == convolution in time domain
amp = ifft(f(:,1).*conj(f(:,2)));

%=== extract portion of the FFT relating to the time lags of interest
lag = 1:N_LAG;
amp = real(amp(lag))/(autocorr_ecg*autocorr_abp);

end