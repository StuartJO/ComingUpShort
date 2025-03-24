
load('Hansen_networks.mat')
A = adj{1};

RewireOrder = {'random','shortest','longest','shortest','longest'};
DistMatch = {'none','prob','prob','invprob','invprob'};
%colors = {'g','r','b','k','c'};

colors = [0, 20, 39; 112, 141, 129; 244, 213, 141;191, 6, 3;141, 8, 1]./255;
rewiring_name = {'Random','Short-similar','Long-similar','Short-dissimilar','Long-dissimilar'};

KScolors = lines(4);
load('Hansen400_rewirings.mat')
GlobalA = [efficiency_bin(A) diffusion_efficiency(A) modularity_Q(A) transitivity_bu(A)];

FtypeName = {'max(KS)','max(RMSE)','max(1-r)','TND','TFdiff'};
FtypeName = {'\it{max(KS)}','\it{max(RMSE)}','\it{max(r_d )}','\it{TND}','\it{TF_{diff}}'};
EdgeOvlp = mean(EdgeOverlap{1})';

TopoType = cell(1,5);
for j = 1:5
TopoType{j} = zeros(1973,5);
TopoType{j}(:,2) = mean(max(TopoMeasures{j}.RMSE(:,:,1:4),[],3));
TopoType{j}(:,5) = mean(TopoMeasures{j}.TF);
GlobalB = TopoMeasures{j}.GlobalTopo;
GolBSz = size(GlobalB);
% reshape(GlobalB,[prod(GolBSz(1:2)) 4]);
TopoDist_ = pdist2(GlobalA,reshape(GlobalB,[prod(GolBSz(1:2)) 4]));
TopoType{j}(:,4) = mean(reshape(TopoDist_,GolBSz(1:2)));
TopoType{j}(:,1) = mean(maxKS{j});
TopoType{j}(:,3) = mean(max(1-TopoMeasures{j}.TopogCorr,[],3));
end

Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);
%ax1{1} = subplot(2,5,[1:3 6:8]);
ax1{1} = subplot(4,4,[1:12]);
for j = 1:5
    plot(1-EdgeOvlp,TopoType{j}(:,1),'Color',colors(j,:),'LineWidth',2)
    hold on
end
%set(gca,'XDir','reverse')
xlabel('Proportion of connections rewired')
ylabel('\it{max(KS)}')
set(gca,'FontSize',18)

lgd = legend(rewiring_name,'Location','southoutside','Orientation','horizontal');

subplotInds = [4 5 9 10];
subplotInds = [13:16];
kstype = {'k','c','b','e'};
for i = 1:4
%ax1{1+i} = subplot(2,5,subplotInds(i));
ax1{1+i} = subplot(4,4,subplotInds(i));
    for j = 1:5
        plot(1-EdgeOvlp,mean(squeeze(KS{j}(:,:,i))),'Color',colors(j,:),'LineWidth',2)
        hold on
    end
    ylim([0 1])
    %set(gca,'XDir','reverse')

xlabel('Proportion of connections rewired')
ylabel(['\it{KS_',kstype{i},'}'])
set(gca,'FontSize',12)
end

AddLetters2Plots(ax1, 'HShift', -0.05, 'VShift', -0.05,'FontSize',24)

exportgraphics(gcf,'./Fig5.png','Resolution',300)

%AddLetters2Plots(ax1, 'HShift', -0.07, 'VShift', -0.04,'FontSize',24)
%%

Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);
for i = 1:5
    ax{i} = subplot(2,5,i);
    for j = 1:5
        plot(1-EdgeOvlp,TopoType{j}(:,i),'Color',colors(j,:),'LineWidth',2)
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
    
    ax{i+5} = subplot(2,5,i+5);
    for j = 1:5
        plot(1-EdgeOvlp,TopoType{i}(:,j),'Color',lines_cmap(j,:),'LineWidth',2)
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
AddLetters2Plots(ax, 'HShift', -0.03, 'VShift', -0.04)

exportgraphics(gcf,'./RewiringAltMeasures.png','Resolution',300)

%%

%EvalMeasures_ = EvalMeasures([1 4 5 3 6]);

Letters={'A','B','C','D','E'};
lines_cmap = lines(5);
figure('Position',[1	49	1382.40000000000	747.200000000000]);
for i = 1:5
    ax{i} = subplot(2,5,i);
    for j = 1:5
        plot(1-EdgeOvlp,TopoType{j}(:,i)./max(TopoType{j}(end,i)),'Color',colors(j,:),'LineWidth',2)
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
    
    ax{i+5} = subplot(2,5,i+5);
    for j = 1:5
        plot(1-EdgeOvlp,TopoType{i}(:,j)./max(TopoType{i}(end,j)),'Color',lines_cmap(j,:),'LineWidth',2)
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
AddLetters2Plots(ax, 'HShift', -0.03, 'VShift', -0.04)

exportgraphics(gcf,'./RewiringAltMeasuresNorm.png','Resolution',300)