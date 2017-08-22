function [ h ] = plotann( data, header, fs, varargin )
%PLOTANN	Plot a subsection of a signal with anns
% [ h ] = plotann( data, header, fs, ann, sqi )
% 
% ann - cell array, annotations in seconds
% sqi - cell array, one annotation per second

%	Copyright 2014 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 25-Aug-2014 13:43:30
%	Contact: alistairewj@gmail.com

%% plotting parameters
col = [0.904705882352941,0.191764705882353,0.198823529411765;0.294117647058824,0.544705882352941,0.749411764705882;0.371764705882353,0.717647058823529,0.361176470588235;1,0.548235294117647,0.100000000000000;0.955000000000000,0.894647058823529,0.472176470588235;0.685882352941177,0.403529411764706,0.241176470588235;0.971764705882353,0.555294117647059,0.774117647058824];
col = [col; col*0.9];
marker = {'o','o','x','x','+','+'};
ms = 16;
%% Define the segment of interest
% T_START = 210*360+1; T_END = min(T_START + 30*360,size(data,1));
T_START = 1; T_END = size(data,1);

%% load in data
[ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header);

if ~isempty(idxECG)
    ecg = data(T_START:T_END,idxECG(1));
else
    ecg = [];
end

if ~isempty(idxABP)
    sig2 = 'Blood pressure';
    bp = data(T_START:T_END,idxABP(1));
elseif ~isempty(idxPPG)
    sig2 = 'Photoplethysmogram';
    bp = data(T_START:T_END,idxPPG(1));
    
elseif ~isempty(idxSV)
    sig2 = 'Stroke volume';
    bp = data(T_START:T_END,idxSV(1));
else 
    bp = [];
    sig2 = [];
end
t = (1:numel(ecg))';
fs = fs(1);

if nargin > 4
   ann = varargin{1};
   sqi = varargin{2};
elseif nargin > 3
    ann = varargin{1};
    sqi = {};
end
%% subsample ANNs
if ~isempty(ann)
    if isnumeric(ann)
        ann = {ann};
    end
    
    for f=1:numel(ann)
        ann{f} = ceil(ann{f}*fs);
        ann{f} = ann{f}(ann{f} >= T_START & ann{f} < T_END); 
        ann{f} = ann{f} - T_START;
    end
else
    ann = {};
end

%% subsample SQIs

if ~isempty(sqi)
    if isnumeric(sqi)
        sqi = {sqi};
    end
    sqi_time = floor(T_START/fs)+1:floor(T_END/fs);
    sqi_time = sqi_time - floor(T_START/fs);
    
    for f=1:numel(sqi)
        sqi{f} = sqi{f}(floor(T_START/fs)+1:floor(T_END/fs)); 
    end
else
    sqi_time = [];
    sqi = {};
end

%=== rescale ecg+bp for plotting
ecg = ecg - min(ecg);
ecg = ecg ./ max(ecg);
ecg = ecg*0.95 + 1.05;

bp = bp - min(bp);
bp = bp ./ max(bp);
bp = bp*0.95 + 3.05;
t = t./max(fs);

%% plot final figure
try
    set(0,'CurrentFigure',1);
catch
    figure(1);
end

clf; hold all;
if ~isempty(ecg)
    plot(t,ecg,'-','linewidth',2,'color',col(1,:));
else
    % dummy plot
    plot(t(1),0,'-','linewidth',2,'color',col(1,:));
end
if ~isempty(bp)
    plot(t,bp,'-','linewidth',2,'color',col(2,:));
else
    plot(t(1),0,'-','linewidth',2,'color',col(2,:));
end

% %=== plot patches for annotation locations
% for f=1:numel(ann.atr)
%     curr_ann = ann.atr(f);
%     patch([-0.15,-0.15,0.15,0.15]+t(curr_ann), [1.05,4.5,4.5,1.05], 0.6*ones(1,3),...
%         'edgecolor','none','facealpha',0.5,'handlevisibility','off');
% end

% %=== plot patches for annotation locations
% for f=1:numel(ann.gqrs)
%     curr_ann = ann.gqrs(f);
%     patch([-0.05,-0.05,0.05,0.05]+t(curr_ann), [1.05,4.5,4.5,1.05], 0.1*ones(1,3),...
%         'edgecolor','none','facealpha',0.5,'handlevisibility','off');
% end

%=== plot ann on ECG
fn_plot = strcat(repmat({'ann'},1,numel(ann)),arrayfun(@num2str,1:numel(ann),'UniformOutput',false));
for f=1:numel(fn_plot)
    curr_ann = ann{f};
    % reformat the marker for clarity
    if any(strcmp({'o','^'},marker{f})); 
        mfill = col(f+2,:); ms = 12; 
    else
        mfill='none'; ms = 16; 
    end
    
    switch fn_plot{f}
        % plot on the abp
        case 'wabp'
            plot(t(curr_ann),2*ones(numel(curr_ann),1) + f / (numel(fn_plot)+1),marker{f},...
                'linewidth',2,'markersize',ms,'color',col(f+2,:),...
                'markerfacecolor',mfill);
        otherwise
            % plot on the ECG
            plot(t(curr_ann),2*ones(numel(curr_ann),1) + f / (numel(fn_plot)+1),marker{f},...
                'linewidth',2,'markersize',ms,'color',col(f+2,:),...
                'markerfacecolor',mfill);
    end
end

% plot SQI
if ~isempty(sqi) && numel(sqi)>0
    for f = 1:numel(sqi)
        % reformat the marker for clarity
        if any(strcmp({'o','^'},marker{f})); 
            mfill = col(f,:); ms = 12; 
        else
            mfill='none'; ms = 16; 
        end
        
        plot(sqi_time, sqi{f},marker{f},...
                'linewidth',2,'markersize',ms,'color',col(f,:),...
                'markerfacecolor',mfill);
    end
end

set(gca,'YLim',[-0.1,4.5]);
% set(gca,'XLim',[0,11]);
hleg=legend(['ECG','BP',fn_plot]);
set(hleg,'Box','Off');

% add some labels to the y-axis for easy reading
set(gca,'YTick',[0:0.5:1,1.5,2 + ((1:numel(fn_plot))/(numel(fn_plot)+1)),3.8],...
    'YTickLabel',['0%','SQI','100%','ECG',fn_plot,sig2]);
xlabel('Time (seconds)');

% plot ground truth if available
% plot(t(ann.atr),ecg(ann.atr),'+','linewidth',2,'markersize',ms,'color',[0,0,0]);

grid on;
end

