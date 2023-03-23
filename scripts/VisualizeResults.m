%% I NEVER WANNA WORK BUT ALSO LIKE BUY ME COOL BIKES
%
%
%#ok<*SAGROW> 

close all force
clear
clc

load('perfeval_ec.mat')
target = [labels{:}];

%% roc and freeman (or is it really southcott?)
%
% loop over networks
tic
for ii = 1:size(NET_OUT,2)
    disp('yoooo look at these results')

    test = [NET_OUT{:,ii}];
    score = [NET_SCR{:,ii}];

    % roc
    R = roc(target,test);
    TPR(ii) = R.TPR;
    FNR(ii) = R.FNR;
    ACC(ii) = R.ACC;
    HFA(ii) = R.HFA;

    % freeman
    F = dynamicMeasures(target,test);
    FEC{ii} = F.FEC;
    OVER{ii} = F.OVER;
    MSC(ii) = F.MSC;
    NDS(ii) = F.NDS;


    % auc
    [X,Y,~,A] = perfcurve(target,score,1);
    AUC1(ii) = A;
    xFPR1{ii} = X;
    yTPR1{ii} = Y;

end

% loop over vads
for ii = 1:size(VAD_OUT,2)
    disp('yoooo look at these VADs')

    test = [VAD_OUT{:,ii}];

    % roc
    R = roc(target,test);
    TPRd(ii) = R.TPR;
    FNRd(ii) = R.FNR;
    ACCd(ii) = R.ACC;
    HFAd(ii) = R.HFA;
    HFAd(ii) = R.HFA;

    % freeman
    F = dynamicMeasures(target,test);
    FECd{ii} = F.FEC;
    OVERd{ii} = F.OVER;
    MSCd(ii) = F.MSC;
    NDSd(ii) = F.NDS;

end
toc

%% gatherup
% 
TP1 = [TPR TPRd] * 100;
FN1 = [FNR FNRd] * 100;
AC1 = [ACC ACCd] * 100;
HF1 = [HFA HFAd] * 100;

FC1 = mean(cell2mat([FEC FECd]),1,'omitnan');
OV1 = mean(cell2mat([OVER OVERd]),1,'omitnan');
MC1 = [MSC MSCd];
NS1 = [NDS NDSd];
 
%%
clearvars -except TP1 FN1 AC1 HF1 FC1 OV1 MC1 NS1 xFPR1 yTPR1 AUC1

load('perfeval_wsa.mat')
target = [labels{:}];

%% roc and freeman (or is it really southcott?)
%
% loop over networks
tic
for ii = 1:size(NET_OUT,2)
    disp('yoooo look at these results')

    test = [NET_OUT{:,ii}];
    score = [NET_SCR{:,ii}];

    % roc
    R = roc(target,test);
    TPR(ii) = R.TPR;
    FNR(ii) = R.FNR;
    ACC(ii) = R.ACC;
    HFA(ii) = R.HFA;

    % freeman
    F = dynamicMeasures(target,test);
    FEC{ii} = F.FEC;
    OVER{ii} = F.OVER;
    MSC(ii) = F.MSC;
    NDS(ii) = F.NDS;

    
    % auc
    [X,Y,~,A] = perfcurve(target,score,1);
    AUC2(ii) = A;
    xFPR2{ii} = X;
    yTPR2{ii} = Y;

end

% loop over vads
for ii = 1:size(VAD_OUT,2)
    disp('yoooo look at these VADs')

    test = [VAD_OUT{:,ii}];

    % roc
    R = roc(target,test);
    TPRd(ii) = R.TPR;
    FNRd(ii) = R.FNR;
    ACCd(ii) = R.ACC;
    HFAd(ii) = R.HFA;
    HFAd(ii) = R.HFA;

    % freeman
    F = dynamicMeasures(target,test);
    FECd{ii} = F.FEC;
    OVERd{ii} = F.OVER;
    MSCd(ii) = F.MSC;
    NDSd(ii) = F.NDS;

end
toc

%% gatherup
% 
TP2 = [TPR TPRd] * 100;
FN2 = [FNR FNRd] * 100;
AC2 = [ACC ACCd] * 100;
HF2 = [HFA HFAd] * 100;

FC2 = mean(cell2mat([FEC FECd]),1,'omitnan');
OV2 = mean(cell2mat([OVER OVERd]),1,'omitnan');
MC2 = [MSC MSCd];
NS2 = [NDS NDSd];

%%
picpath = '/Users/akittel/Documents/dutten/yr7 thesis/latex/pics/';
W = 18.2;
L = 10;

% Colormap
n = 5;
cmap = colormapVoicebox(n);

% Marker
strMarker = {'^' 'v' '<' '>' 'o' 's' 'd' 'p' 'h' '*' '+' 'x' '^' 'v' '<' '>' 'o' 's' 'd' 'p' 'h' '*' '+' 'x'};

% Legend
netIds = {'biLSTM MFCC','biLSTM mel-spectr.','LSTM MFCC'}.';
map = [netIds;{'Ghosh LTSV';'Graf COMB'}];

%% 
% plot1
% 
figure(1)
t=tiledlayout(3,2,'TileSpacing','compact');

% roc
ax = nexttile(1);
vals = [TP2;FN2;HF2;AC2];
xTickNames = ["TPR","FNR","HFA","ACC"];

b = bar(1:4,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('%-score')
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
set(ax,'YTick',[0 25 50 75 100])
ylim([0 100])
title(ax,'WSA test data')
text(0.15,0.88,'(a)','Units','normalized','FontSize',8)

% freeman NDS and MSC
ax = nexttile(3);
vals = [MC2;NS2];
xTickNames = ["MSC","NDS"];

b = bar(1:2,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('%-score')
yticks([0 5 10 15])
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
ylim([0 17])
text(0.02,0.88,'(c)','Units','normalized','FontSize',8)


% freeman FEC and OVER
ax = nexttile(5);
vals = [FC2;OV2];
xTickNames = ["FEC","OVER"];

b = bar(1:2,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('Avg. duration (ms)')
yticks([0 1 2 3 4])
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
text(0.02,0.88,'(e)','Units','normalized','FontSize',8)


% roc
ax = nexttile(2);
vals = [TP1;FN1;HF1;AC1];
xTickNames = ["TPR","FNR","HFA","ACC"];

b = bar(1:4,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('%-score')
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
set(ax,'YTick',[0 25 50 75 100])
ylim([0 100])
title(ax,'EC dataset')
text(0.02,0.88,'(b)','Units','normalized','FontSize',8)


% freeman NDS and MSC
ax = nexttile(4);
vals = [MC1;NS1];
xTickNames = ["MSC","NDS"];

b = bar(1:2,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('%-score')
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
ylim([0 30])
text(0.02,0.88,'(d)','Units','normalized','FontSize',8)


% freeman FEC and OVER
ax = nexttile(6);
vals = [FC1;OV1];
xTickNames = ["FEC","OVER"];

b = bar(1:2,vals');
set(b, {'DisplayName'}, map);
colororder(cmap)
grid on; box on;
ylabel('Avg. duration (ms)')
yticks([0 2.5 5 7.5 10])
set(ax,'XTickLabels',xTickNames,'FontSize', 8)
text(0.02,0.88,'(f)','Units','normalized','FontSize',8)


lg = legend(ax);
lg.Interpreter = 'none';
lg.Layout.Tile = 'South';
lg.NumColumns = 2;


set(gcf, ...
    'Units', 'centimeter');
set(gcf, ...
    'Position', [30 5 W L]);
exportgraphics(gcf,[picpath,'results.eps'],'BackgroundColor','none','ContentType','vector')

% plot2
% 
% strLegend = cell(numel(netIds),1);
% figure(2)
% tiledlayout(2,1)
% nexttile
% hold on
% grid on; box on;
% for ii = 1:size(NET_OUT,2)
%     h = plot(xFPR2{1,ii},yTPR2{1,ii});
%     set(h,'LineWidth',1,'color',cmap(ii,:));
% end
% plot(0:0.1:1,0:0.1:1,'--','Color','#5A5A5A')
% title('WSA test data')
% ylabel('TPR')
% xlabel('FPR')
% set(gca,'FontSize', 8)
% 
% hold off
% 
% nexttile
% hold on
% grid on; box on;
% for ii = 1:size(NET_OUT,2)
%     h = plot(xFPR1{1,ii},yTPR1{1,ii});
%     set(h,'LineWidth',1,'color',cmap(ii,:));
%     strLegend{ii} = netIds{ii};
% end
% set(gca,'FontSize', 8)
% plot(0:0.1:1,0:0.1:1,'--','Color','#5A5A5A')
% title('EC dataset')
% lg = legend(strLegend);
% ylabel('TPR')
% xlabel('FPR')
% hold off
% 
% % lg = legend(ax);
% lg.Interpreter = 'none';
% lg.Layout.Tile = 'South';
% 
% set(gcf, ...
%     'Units', 'centimeter');
% set(gcf, ...
%     'Position', [15 5 9 8]);
% exportgraphics(gcf,[picpath,'auc.eps'],'BackgroundColor','none','ContentType','vector')