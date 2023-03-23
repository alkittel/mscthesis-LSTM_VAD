%#ok<*SAGROW>
%#ok<*PFBNS>

close all force
clear
clc

load('/Users/akittel/Documents/dutten/yr7 thesis/src/MAT-files/results/final/perfeval_snr.mat')

nConditions = size(conditions,1);
nNets = size(NET_OUT,2);
nVads = size(VAD_OUT,2);
n = nNets + nVads;

idxMixtures = conditions(:,1);
idxSNR      = conditions(:,2);
idxNoise    = conditions(:,3);
snrdBVec    = snrdB(idxSNR);

%% roc
%
tic
% loop over conditions
parfor ii = 1:nConditions
    target = labels{ii};
    % loop over networks
    for jj = 1:nNets
        disp('yoooo look at these results')

        test = NET_OUT{ii,jj};
        score = NET_SCR{ii,jj};

        % roc
        R = roc(target,test);
        TPR(ii,jj) = R.TPR;
        FNR(ii,jj) = R.FNR;
        ACC(ii,jj) = R.ACC;
        HFA(ii,jj) = R.HFA;

        % freeman
        F = dynamicMeasures(target,test);
        FEC{ii,jj} = F.FEC;
        OVER{ii,jj} = F.OVER;
        MSC(ii,jj) = F.MSC;
        NDS(ii,jj) = F.NDS;

        % auc
        [~,~,~,A] = perfcurve(target,score,1);
        AUC(ii,jj) = A;
    end

    % loop over vads
    for kk = 1:nVads
        disp('yoooo look at these VADs')

        test = VAD_OUT{ii,kk};

        % roc
        R = roc(target,test);
        TPRd(ii,kk) = R.TPR;
        FNRd(ii,kk) = R.FNR;
        ACCd(ii,kk) = R.ACC;
        HFAd(ii,kk) = R.HFA;

        % freeman
        F = dynamicMeasures(target,test);
        FECd{ii,kk} = F.FEC;
        OVERd{ii,kk} = F.OVER;
        MSCd(ii,kk) = F.MSC;
        NDSd(ii,kk) = F.NDS;
    end
end
toc

%% snr analysis (see experiment_vad)
%
idxSNR = conditions(:,2);
nSNRs = idxSNR(end);

% prealloc
TP = zeros(nSNRs,n);
FN = zeros(nSNRs,n);
AC = zeros(nSNRs,n);
HF = zeros(nSNRs,n);

MC = zeros(nSNRs,n);
NS = zeros(nSNRs,n);
FC = zeros(nSNRs,n);
OV = zeros(nSNRs,n);

FEC2 = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), FEC, 'uni', false));
FEC2d = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), FECd, 'uni', false));

OVER2 = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), OVER, 'uni', false));
OVER2d = cell2mat(cellfun(@(x) mean(x,1,'omitnan'), OVERd, 'uni', false));

AU = zeros(nSNRs,nNets);


for ii = 1:nSNRs
    TP(ii,:) = [100 * mean(TPR(idxSNR == ii,:),'omitnan') ...
        100 * mean(TPRd(idxSNR == ii,:),'omitnan')];

    FN(ii,:) = [100 * mean(FNR(idxSNR == ii,:),'omitnan') ...
        100 * mean(FNRd(idxSNR == ii,:),'omitnan')];

    AC(ii,:) = [100 * mean(ACC(idxSNR == ii,:),'omitnan') ...
        100 * mean(ACCd(idxSNR == ii,:),'omitnan')];

    HF(ii,:) = [100 * mean(HFA(idxSNR == ii,:),'omitnan') ...
        100 * mean(HFAd(idxSNR == ii,:),'omitnan')];

    MC(ii,:) = [mean(MSC(idxSNR == ii,:),'omitnan') ...
        mean(MSCd(idxSNR == ii,:),'omitnan')];

    NS(ii,:) = [mean(NDS(idxSNR == ii,:),'omitnan') ...
        mean(NDSd(idxSNR == ii,:),'omitnan')];

    FC(ii,:) = [mean(FEC2(idxSNR == ii,:),'omitnan') ...
        mean(FEC2d(idxSNR == ii,:),'omitnan')];

    OV(ii,:) = [mean(OVER2(idxSNR == ii,:),'omitnan') ...
        mean(OVER2d(idxSNR == ii,:),'omitnan')];

    AU(ii,:) = mean(AUC(idxSNR == ii,:),'omitnan');
end

% Colormap
cmap = colormapVoicebox(n);

% Marker
strMarker = {'^' 'v' '<' '>' 'o' 's' 'd' 'p' 'h' '*' '+' 'x' '^' 'v' '<' '>' 'o' 's' 'd' 'p' 'h' '*' '+' 'x'};

% Legend
netIds = {'biLSTM MFCC','biLSTM mel-spectr.','LSTM MFCC'}.';
map = [netIds;{'Ghosh LTSV';'Graf COMB'}];
strLegend = cell(n,1);

%% plots
picpath = '/Users/akittel/Documents/dutten/yr7 thesis/latex/pics/';
W = 20;
L = 8;

close all

figure(1);
tiledlayout(4,2,'TileSpacing','compact')
nexttile(1)
hold on;
for ii = 1 : n
    h = plot(snrdB,TP(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
    strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('TPR (%)')
xlim([-11 16])
set(gca,'FontSize', 8)

nexttile(2)
hold on;
for ii = 1 : n
    h = plot(snrdB,FN(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
%     strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('FNR (%)')
xlim([-11 16])
set(gca,'FontSize', 8)

nexttile(3)
hold on;
for ii = 1 : n
    h = plot(snrdB,HF(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
%     strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('HFA (%)')
xlim([-11 16])
set(gca,'FontSize', 8)

nexttile(4)
hold on;
for ii = 1 : n
    h = plot(snrdB,AC(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
%     strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('ACC (%)')
xlim([-11 16])
set(gca,'FontSize', 8)

nexttile(5)
hold on;
for ii = 1 : n
    h = plot(snrdB,MC(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
    strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('MSC (%)')
xlim([-11 16])
ylim([0 15])
yticks([0 5 10 15])
% ylim([0 11])
set(gca,'FontSize', 8)

nexttile(6)
hold on;
for ii = 1 : n
    h = plot(snrdB,NS(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
    strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('NDS (%)')
xlim([-11 16])
ylim([0 8])
yticks([0 2 4 6 8])
% ylim([0 11])
set(gca,'FontSize', 8)

nexttile(7)
hold on;
for ii = 1 : n
    h = plot(snrdB,FC(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
    strLegend{ii} = map{ii};
end
grid on; box on;
ylabel('FEC (ms)')
xlabel('SNR (dB)')
ylim([0 2.5])
yticks([0 0.5 1 1.5 2 2.5])
xlim([-11 16])
% ylim([0 3.2])
set(gca,'FontSize', 8)

nexttile(8)
hold on;
for ii = 1 : n
    h = plot(snrdB,OV(:,ii));
    set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
        'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
    strLegend{ii} = map{ii};
end
grid on; box on;
xlabel('SNR (dB)')
ylabel('OVER (ms)')
ylim([0 2.5])
yticks([0 0.5 1 1.5 2 2.5])
xlim([-11 16])
set(gca,'FontSize', 8)

l = legend(strLegend);
l.Layout.Tile = 'South';
l.NumColumns = 3;

set(gcf, ...
    'Units', 'centimeter');
set(gcf, ...
    'Position', [22 15 W 13]);
exportgraphics(gcf,[picpath,'snr_rocdyn.eps'],'BackgroundColor','none','ContentType','vector')


% strLegend = cell(nNets,1);
% 
% figure(2);
% tiledlayout(1,1)
% nexttile
% hold on;
% for ii = 1 : nNets
%     h = plot(snrdB,AU(:,ii));
%     set(h,'LineWidth',1,'Marker',strMarker{ii},'MarkerSize',6,...
%         'color',cmap(ii,:),'markerfacecolor',cmap(ii,:));
%      strLegend{ii} = map{ii};
% end
% l = legend(strLegend);
% l.Location = 'SouthOutside';
% l.NumColumns = 1;
% grid on; box on;
% xlabel('SNR (dB)')
% ylabel('AUC')
% xlim([-11 16])
% set(gca,'FontSize', 8)
% 
% set(gcf, ...
%     'Units', 'centimeter');
% set(gcf, ...
%     'Position', [2 15 9 5]);
% exportgraphics(gcf,[picpath,'snr_auc.eps'],'BackgroundColor','none','ContentType','vector')
