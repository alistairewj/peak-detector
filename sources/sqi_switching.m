function [qrs_final, leadUsed] = sqi_switching(qrs, sqi, tsqi, THR, delay, plotFlag)
%SQI_SWITCHING	Switch QRS detections based upon SQIs
%	[ qrs_final, leadUsed] ] = SQI_SWITCHING(qrs,sqi,tsqi) switches between 
%	multiple inputs by comparing a signal quality index (SQI) evaluated on 
%	each signal. That is, if there are two signals, the function will pick:
%		i) the first signal with SQI > 0.8
%		ii) the signal with the highest SQI
%
%	[ qrs_final, leadUsed] ] = SQI_SWITCHING(qrs,sqi,tsqi,THR) allows
%	specifying the threshold at which to start switching (above 0.8).
%	The first signal with an SQI above this threshold is used.
%
%	Inputs:
%   	qrs 		-	1xD cell of QRS detections (in seconds)
%   	sqi 		- 	1xD cell of SQI values or WxD matrix. 
%				Note: SQIs for signals must be measured at the same time.
%				Otherwise this function will not work.
%   	
%		tsqi 		- Wx1 matrix of starting times for each SQI value.
%		THR 		- Scalar representing the SQI switching threshold
%		plotFlag 	- Flag for debugging plot
%
%	Outputs:
%		qrs_final	- Final QRS output after SQI switching
%		leadUsed	- Wx1 matrix indicating the lead used in the switch
%
%	See also CHALLENGE ECGSQI

%	Copyright 2014 Alistair Johnson

%=== Preamble - check inputs and get constants
if nargin<4
    THR = 0.8;
end
if nargin<5
    delay = 0.2; % delay between ECG and ABP, also known as pulse transit time
end
if nargin<6
	plotFlag = 0;
end

%=== the SQIs should be the same size - sometimes input as a cell though
if iscell(sqi)
    sqi = horzcat(sqi{:});
end
[W,D] = size(sqi); % number of windows // number of signals

%=== 1) Loop through the SQIs and decide which lead to use
leadUsed = zeros(W,1);
for w=1:W
    %=== if all the leads have reduced signal quality
    if all(sqi(w,:) <= THR)
        [~,leadUsed(w)] = max(sqi(w,:)); % we pick the best one
    else
        % we pick the first with SQI > THR
        leadUsed(w) = find(sqi(w,:)>THR,1);
    end
end
%=== 2) Loop through the signal and allocate the final peaks
qrs_final = cell(W,1);

% the first SQI should always start from the start of the signal
w=1;
while w <= W
    idxWin = tsqi(w);
    if w==W
        idxWinNext = Inf; % take everything until end of signal
    else
        idxWinNext = tsqi(w+1);
    end
    idxQRS = qrs{leadUsed(w)} > idxWin & qrs{leadUsed(w)} <= idxWinNext;
    qrs_final{w} = qrs{leadUsed(w)}(idxQRS);
    
    %=== remove possible double detections due to switching
    if w > 1 && leadUsed(w)~=leadUsed(w-1) && ~isempty(qrs_final{w}) && ~isempty(qrs_final{w-1})
        % can happen on the edge of windows
        % 100ms window for merging
        if qrs_final{w-1}(end) > (qrs_final{w}(1)-0.100)
            % average the two locations
            qrs_final{w-1}(end) = (qrs_final{w-1}(end) + qrs_final{w}(1))/2;
            % delete the other beat
            qrs_final{w}(1) = [];
        end
    end
    
    %=== add in possible missed detections due to switching
    if w > 1 && leadUsed(w)~=leadUsed(w-1)
        % can happen on the edge of windows
        % check for a beat for the other lead 
        beatJustBefore = qrs{leadUsed(w)}(find(qrs{leadUsed(w)}<=idxWin,1,'last'));
        beatJustAfter = qrs{leadUsed(w-1)}(find(qrs{leadUsed(w-1)}>idxWin,1,'first'));
        
        if ~isempty(beatJustBefore) && ~isempty(beatJustAfter) && ...
                beatJustBefore >= (tsqi(w)-delay(leadUsed(w))) && beatJustAfter < (tsqi(w)+delay(leadUsed(w)))
            qrs_final{w} = [beatJustBefore; qrs_final{w}];
        end
    end
    w=w+1;
end
qrs_final = vertcat(qrs_final{:});
qrs_final = qrs_final(:);

if plotFlag
   col = [0.9047    0.1918    0.1988
    0.2941    0.5447    0.7494
    0.3718    0.7176    0.3612
    1.0000    0.5482    0.1000
    0.9550    0.8946    0.4722
    0.6859    0.4035    0.2412
    0.9718    0.5553    0.7741];
    lineStyle={'-','--',':','-','--',':','-'};
   lw = 3;
   
   figure(1); clf; hold all;
   %=== plot heart rates
   for ss=1:D
        plot(qrs{ss}(1:end-1),60./diff(qrs{ss}),lineStyle{ss},'Color',col(ss,:),'LineWidth',lw);
   end
   %=== plot final heart rate
   plot(qrs_final(1:end-1),60./diff(qrs_final),'k','LineWidth',1);
   
   %=== plot SQI
   for ss=1:D
        plot(tsqi,sqi(:,ss)*100,':','Color',col(ss,:),'LineWidth',lw);
   end
   legendStr = [strcat('HR',arrayfun(@num2str, 1:D, 'UniformOutput',false)),...
       'merged',...
       strcat('SQI',arrayfun(@num2str, 1:D, 'UniformOutput',false))];
   legend(legendStr);
   xlabel('Time (seconds)','FontSize',16);
   set(gca,'FontSize',16);
end

end
















