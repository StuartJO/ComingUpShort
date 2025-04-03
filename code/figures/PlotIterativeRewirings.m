function PlotIterativeRewirings

load('Hansen_networks.mat')
A = adj{1};

RewireOrder = {'random','shortest','longest','shortest','longest'};
DistMatch = {'none','prob','prob','invprob','invprob'};

Fullcolors = [0, 20, 39; 112, 141, 129; 244, 213, 141;191, 6, 3;141, 8, 1]./255;

colors = Fullcolors([1 3 5],:);


rewiring_name = {'Random','Short','Long','Short-dissimilar','Long-dissimilar'};

KScolors = lines(4);

load('Hansen400_rewirings.mat')


FtypeName = {'max\it{(KS)}','max\it{(RMSE)}','max\it{(r_d )}','\it{TND}','\it{TF_{diff}}'};
EdgeOvlp = mean(EdgeOverlap{1})';

FitStats = cell(1,3);
GlobalB = cell(1,3);
for j = 1:3
FitStats{j} = zeros(1973,5);
FitStats{j}(:,2) = mean(FitMeasures{j}.maxRMSE);
FitStats{j}(:,5) = mean(FitMeasures{j}.TF);
FitStats{j}(:,4) = mean(FitMeasures{j}.TND);
FitStats{j}(:,1) = mean(FitMeasures{j}.maxKS);
FitStats{j}(:,3) = mean(FitMeasures{j}.maxRd);
end

Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);

ax1{1} = subplot(4,4,1:12);

hold on

load('BestMdls_GNM_maxKS_Add_exponential.mat', 'maxKS')

BestMdlmaxKS = min(mean(maxKS,2));

plot([0 1],[BestMdlmaxKS BestMdlmaxKS],'k','LineStyle','-')

for j = 1:3
    plot(1-EdgeOvlp,FitStats{j}(:,1),'Color',colors(j,:),'LineWidth',2)
    hold on
end

xlabel('Proportion of connections rewired')
ylabel('max\it{(KS)}')
set(gca,'FontSize',18)

lgd = legend([{'Best GNM fit'},rewiring_name],'Location','southoutside','Orientation','horizontal');

subplotInds = [4 5 9 10];
subplotInds = [13:16];
kstype = {'k','c','b','e'};
for i = 1:4

ax1{1+i} = subplot(4,4,subplotInds(i));
    for j = 1:3
        plot(1-EdgeOvlp,mean(squeeze(FitMeasures{j}.KS(:,:,i))),'Color',colors(j,:),'LineWidth',2)
        hold on
    end
    ylim([0 1])
    %set(gca,'XDir','reverse')

xlabel('Proportion of connections rewired')
ylabel(['\it{KS_',kstype{i},'}'])
set(gca,'FontSize',12)
end

AddLetters2Plots(ax1, 'HShift', -0.05, 'VShift', -0.05,'FontSize',24)


exportgraphics(gcf,'./figures/Figure5.png','Resolution',300)

%AddLetters2Plots(ax1, 'HShift', -0.07, 'VShift', -0.04,'FontSize',24)
%%
Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);

subplotINDStop = [1:3;4:6;7:9;10:12;13:15];
subplotINDSbottom = [16:20; 21:25; 26:30];


subplotINDStop = [1:4; 6:9; 11:14; 16:19; 21:24];
subplotINDSbottom = [25:31; 33:39; 42:48];

%tiledlayout(2,15)

for i = 1:5
    ax{i} = subplot(2,24,subplotINDStop(i,:));
    %ax{i} = nexttile([1 3]);
    for j = 1:3
        plot(1-EdgeOvlp,FitStats{j}(:,i),'Color',colors(j,:),'LineWidth',2)
        hold on
    end
    %set(gca,'XDir','reverse')
    xlabel('Proportion of connections rewired')
    ylabel(FtypeName{i})
    if i == 3
    lgd = legend(rewiring_name,'Location','northoutside','Orientation','horizontal');    
    lgdPosition = lgd.Position;
    lgd.Position(1) = (1-lgdPosition(3))/2;
    axPos = ax{3}.OuterPosition;
    lgd.Position(2) = axPos(2)-(lgdPosition(4)*1.5);
    end
end

for i = 1:3
    ax{i+5} = subplot(2,24,subplotINDSbottom(i,:));
    %ax{i+5} = nexttile([1 5]);
    for j = 1:5
        plot(1-EdgeOvlp,FitStats{i}(:,j),'Color',lines_cmap(j,:),'LineWidth',2)
        hold on
    end
    %set(gca,'XDir','reverse')
    xlabel('Proportion of connections rewired')
    ylabel('Statistic')
    title(rewiring_name{i})
    if i == 3
    lgd2 = legend(FtypeName,'Location','northoutside','Orientation','horizontal');
    lgdPosition = lgd2.Position;
    lgd2.Position(1) = (1-lgdPosition(3))/2;
    axPos = ax{i+5}.OuterPosition;
    lgd2.Position(2) = axPos(2)-(lgdPosition(4)*1.5);
    end
end

Ax8_x = ax{8}.Position(1);
Ax7_x = ax{7}.Position(1);
Ax6_x = ax{6}.Position(1);

ax{7}.Position(1) = (Ax8_x+Ax6_x)/2;

AddLetters2Plots(ax, 'HShift', -0.03, 'VShift', -0.04)

%exportgraphics(gcf,'./RewiringAltMeasures.png','Resolution',300)
exportgraphics(gcf,'./figures/S15.png','Resolution',300)

%%
Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);

subplotINDStop = [1:3;4:6;7:9;10:12;13:15];
subplotINDSbottom = [16:20; 21:25; 26:30];


subplotINDStop = [1:4; 6:9; 11:14; 16:19; 21:24];
subplotINDSbottom = [25:31; 33:39; 42:48];

%tiledlayout(2,15)

for i = 1:5
    ax{i} = subplot(2,24,subplotINDStop(i,:));
    %ax{i} = nexttile([1 3]);
    for j = 1:3
        plot(1-EdgeOvlp,FitStats{j}(:,i)./max(FitStats{j}(end,i)),'Color',colors(j,:),'LineWidth',2)
        hold on
    end
    %set(gca,'XDir','reverse')
    xlabel('Proportion of connections rewired')
    ylabel(FtypeName{i})
    if i == 3
    lgd = legend(rewiring_name,'Location','northoutside','Orientation','horizontal');    
    lgdPosition = lgd.Position;
    lgd.Position(1) = (1-lgdPosition(3))/2;
    axPos = ax{3}.OuterPosition;
    lgd.Position(2) = axPos(2)-(lgdPosition(4)*1.5);
    end
end

for i = 1:3
    ax{i+5} = subplot(2,24,subplotINDSbottom(i,:));
    %ax{i+5} = nexttile([1 5]);
    for j = 1:5
        plot(1-EdgeOvlp,FitStats{i}(:,j)./max(FitStats{i}(end,j)),'Color',lines_cmap(j,:),'LineWidth',2)
        hold on
    end
    %set(gca,'XDir','reverse')
    xlabel('Proportion of connections rewired')
    ylabel('Statistic')
    title(rewiring_name{i})
    if i == 3
    lgd2 = legend(FtypeName,'Location','northoutside','Orientation','horizontal');
    lgdPosition = lgd2.Position;
    lgd2.Position(1) = (1-lgdPosition(3))/2;
    axPos = ax{i+5}.OuterPosition;
    lgd2.Position(2) = axPos(2)-(lgdPosition(4)*1.5);
    end
end

Ax8_x = ax{8}.Position(1);
Ax7_x = ax{7}.Position(1);
Ax6_x = ax{6}.Position(1);

ax{7}.Position(1) = (Ax8_x+Ax6_x)/2;

AddLetters2Plots(ax, 'HShift', -0.03, 'VShift', -0.04)

%exportgraphics(gcf,'./RewiringAltMeasuresNorm.png','Resolution',300)

%%
figure
ax1_{1} = subplot(2,1,1);

for j = 1:3
    plot(1-EdgeOvlp,mean(FitMeasures{j}.DegCorr),'Color',colors(j,:),'LineWidth',2)
    hold on
end
%set(gca,'XDir','reverse')
xlabel('Proportion of connections rewired')
ylabel('Degree correlation')

ax1_{2} = subplot(2,1,2);

for j = 1:3
    %plot(TopoType{j}(:,1),DegCorr{j},'Color',colors(j,:),'LineWidth',2)
    plot(1-EdgeOvlp,FitStats{j}(:,1),'Color',colors(j,:),'LineWidth',2)
    hold on
end
xlabel('Proportion of connections rewired')
ylabel('max\it{(KS)}')
lgd = legend(rewiring_name,'Location','southoutside','Orientation','horizontal');

AddLetters2Plots(ax1_, 'HShift', -0.1, 'VShift', -0.05,'FontSize',18)
%exportgraphics(gcf,'./RewiringDegree.png','Resolution',300)
exportgraphics(gcf,'./figures/S13.png','Resolution',300)