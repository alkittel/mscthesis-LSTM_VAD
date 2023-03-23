%% INITIALIZATION
%#ok<*SAGROW>
%
%

clear
clc

cd('/Users/akittel/Documents/dutten/yr7 thesis')
addpath(genpath(pwd))

T = 300; % WSAMIX
parameters

rng(0)

%% collect validation results
mat1 = matfile('8mfcc_hpval_mix.mat');
mat2 = matfile('16mfcc_hpval_mix.mat');
mat3 = matfile('24mfcc_hpval_mix.mat');

netout1 = mat1.NET_OUT;
netout2 = mat2.NET_OUT;
netout3 = mat3.NET_OUT;
NET_OUT = [netout1 netout2 netout3];

labels = [repmat(mat1.labels,1,10) ...
    repmat(mat2.labels,1,2) ...
    repmat(mat3.labels,1,2)];

dir1 = mat1.dirlist;
dir2 = mat2.dirlist;
dir3 = mat3.dirlist;
dirlist = [dir1;dir2;dir3];
for ii = 1:numel(dirlist)
    netMap{ii} = dirlist(ii).name(1:end-24);
end

%% run
for ii = 1:size(NET_OUT,2)
    disp('yoooo look at these results')
    test = [NET_OUT{:,ii}];
    target = [labels{:,ii}];

    % roc
    R = roc(target,test);
    ACC1(ii) = R.ACC;
    HFA1(ii) = R.HFA;
end


%% INITIALIZATION
%#ok<*SAGROW>
%
%

clearvars -except ACC1 HFA1
clc

cd('/Users/akittel/Documents/dutten/yr7 thesis')
addpath(genpath(pwd))

T = 300; % WSAMIX
parameters

rng(0)

%% collect validation results
mat1 = matfile('8mfcc_hpval_wsa.mat');
mat2 = matfile('16mfcc_hpval_wsa.mat');
mat3 = matfile('24mfcc_hpval_wsa.mat');

netout1 = mat1.NET_OUT;
netout2 = mat2.NET_OUT;
netout3 = mat3.NET_OUT;
NET_OUT = [netout1 netout2 netout3];

labels = [repmat(mat1.labels,1,10) ...
    repmat(mat2.labels,1,2) ...
    repmat(mat3.labels,1,2)];

dir1 = mat1.dirlist;
dir2 = mat2.dirlist;
dir3 = mat3.dirlist;
dirlist = [dir1;dir2;dir3];
for ii = 1:numel(dirlist)
    netMap{ii} = dirlist(ii).name(1:end-24);
end

netMap = cellfun(@(x) strrep(x,'DROPOUT','DP'), netMap,'uni',0);
netMap = cellfun(@(x) strrep(x,'BATCHNORM','BN'), netMap,'uni',0);
netMap = cellfun(@(x) strrep(x,'NONORM','NN'), netMap,'uni',0);


%% run
for ii = 1:size(NET_OUT,2)
    disp('yoooo look at these results')
    test = [NET_OUT{:,ii}];
    target = [labels{:,ii}];

    % roc
    R = roc(target,test);
    ACC2(ii) = R.ACC;
    HFA2(ii) = R.HFA;
end


%% plot
W = 9;

cmap = colormapVoicebox(14);

% roc
figure(1)
tiledlayout(2,1)
vals1 = [HFA1 * 100 ; ACC1 * 100];
vals2 = [HFA2 * 100 ; ACC2 * 100];
xTickNames = ["HFA (%)","ACC (%)"];

ax = nexttile(1);

b = bar(1:2,vals1');
set(b, {'DisplayName'}, netMap.');
colororder(cmap)
grid on
ylim(ax,[60 100])
set(ax,'XTickLabels',xTickNames,'FontSize',8)
title(ax,'Validation metrics, WSA + mixtures')

ax = nexttile(2);

b = bar(1:2,vals2');
set(b, {'DisplayName'}, netMap.');
colororder(cmap)
grid on
ylim(ax,[85 100])
set(ax,'XTickLabels',xTickNames,'FontSize',8)
title(ax,'Validation metrics, WSA only')

lg = legend(ax);
lg.Interpreter = 'none';
lg.Layout.Tile = 'West';
lg.NumColumns = 1;
lg.ItemTokenSize = [15 18];

set(gcf,'renderer','Painters')
set(gcf, ...
    'Units', 'centimeter');
set(gcf, ...
    'Position', [15 5 18.5 6.5]);
exportgraphics(gcf,[picpath,'validation.eps'],'BackgroundColor','none','ContentType','vector')

