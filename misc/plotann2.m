function [ h ] = plotann2( data, header, fs, varargin )
%PLOTANN	Plot a subsection of a signal with anns

%	Copyright 2014 Alistair Johnson

%	$LastChangedBy$
%	$LastChangedDate$
%	$Revision$
%	Originally written on GLNXA64 by Alistair Johnson, 25-Aug-2014 13:43:30
%	Contact: alistairewj@gmail.com

%% plotting parameters
col = [0.904705882352941,0.191764705882353,0.198823529411765;0.294117647058824,0.544705882352941,0.749411764705882;0.371764705882353,0.717647058823529,0.361176470588235;1,0.548235294117647,0.100000000000000;0.955000000000000,0.894647058823529,0.472176470588235;0.685882352941177,0.403529411764706,0.241176470588235;0.971764705882353,0.555294117647059,0.774117647058824];
marker = {'o','x','o','x','+','+'};
ms = 16;
ann_plot_type = 2; % if ==1 then patches are plotted, otherwise symbols.
%% Define the segment of interest
% T_START = 210*360+1; T_END = min(T_START + 30*360,size(data,1));
T_START = 290*fs+1; T_END = 300*fs;
[N,D] = size(data);

%% load in data
[ idxECG, idxABP, idxPPG, idxSV ] = getSignalIndices(header);

%=== create the header
sig_header = [ strcat(repmat({'ECG '}, 1, numel(idxECG)), arrayfun(@num2str, 1:numel(idxECG), 'UniformOutput', false)),...
    strcat(repmat({'ABP '}, 1, numel(idxABP)), arrayfun(@num2str, 1:numel(idxABP), 'UniformOutput', false)),...
    strcat(repmat({'PPG '}, 1, numel(idxPPG)), arrayfun(@num2str, 1:numel(idxPPG), 'UniformOutput', false)),...
    strcat(repmat({'SV '}, 1, numel(idxSV)), arrayfun(@num2str, 1:numel(idxSV), 'UniformOutput', false))];

% subsample + select signals of interest
data = data(T_START:T_END,[idxECG, idxABP, idxPPG, idxSV]);



%% create and subsample annotations
if nargin>2
    ann = varargin;
else
    ann = [];
end

for f=1:numel(ann)
    ann{f} = ann{f}(ann{f} >= T_START & ann{f} < T_END); 
    ann{f} = ann{f} - T_START;
end

%%
%=== rescale data to [0,1] for plotting
data = bsxfun(@minus, data, min(data,[],1));
data = bsxfun(@rdivide, data, max(data,[],1));

idxFix = all(isnan(data),1);
data(:,idxFix) = 0; % if data is all the same number, set it to 0

% Shift the data so it ranges from [1-2], [2-3], etc etc
for d=1:size(data,2)
    data(:,d) = data(:,d) * 0.95 + d + 0.025;
end
%% plot final figure
t = (1:size(data,1))'-1; % time vector
t = t./fs;
h=figure(1); clf; hold all;
for d=1:size(data,2)
    plot(t,data(:,d),'-','linewidth',2,'color',col(d,:));
end

fn_plot = strcat(repmat({'ann'},1,numel(ann)),arrayfun(@num2str,1:numel(ann),'UniformOutput',false));
switch ann_plot_type
    case 1
        % %=== plot patches for annotation locations
        % for f=1:numel(ann.atr)
        %     curr_ann = ann.atr(f);
        %     patch([-0.15,-0.15,0.15,0.15]+t(curr_ann), [1.05,4.5,4.5,1.05], 0.6*ones(1,3),...
        %         'edgecolor','none','facealpha',0.5,'handlevisibility','off');
        % end
        
        %=== plot patches for annotation locations
        for a=1:numel(ann)
            curr_ann_plotted = ann{a};
            for f=1:numel(curr_ann_plotted)
                curr_ann = curr_ann_plotted(f);
                patch([-0.05,-0.05,0.05,0.05]+t(curr_ann), [1.05,D+1,D+1,1.05], col(a,:),...
                    'edgecolor','none','facealpha',0.5,'handlevisibility','off');
            end
        end
        
    case 2
        %=== plot patches for annotation locations
        for a=1:numel(ann)
            plot([t(ann{a}),t(ann{a})]', ([0;D+1]*ones(1,numel(ann{a}))), '--', 'color', col(a,:),...
                'linewidth',1.5,'marker',marker{a},'markersize',16,'handlevisibility','off');
        end
        
    otherwise
        %=== plot ann on ECG
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
                        'markerfacecolor',mfill,'handlevisibility','off');
                otherwise
                    % plot on the ECG
                    plot(t(curr_ann),2*ones(numel(curr_ann),1) + f / (numel(fn_plot)+1),marker{f},...
                        'linewidth',2,'markersize',ms,'color',col(f+2,:),...
                        'markerfacecolor',mfill,'handlevisibility','off');
            end
        end
end
        
%=== dummy plotted signals for the legend call
for a=1:numel(ann)
    plot(t(1), -0.5, '-','color',col(a,:),'linewidth',3,'marker','none');
end

set(gca,'YLim',[-0.1,size(data,2)+1]);
% set(gca,'XLim',[0,11]);
hleg=legend([sig_header,fn_plot]);
set(hleg,'Box','Off');

% add some labels to the y-axis for easy reading
set(gca,'YTick',(1:size(data,2))+0.5,...
    'YTickLabel',sig_header);
xlabel('Time (seconds)');

% plot ground truth if available
% plot(t(ann.atr),ecg(ann.atr),'+','linewidth',2,'markersize',ms,'color',[0,0,0]);

grid on;
end

