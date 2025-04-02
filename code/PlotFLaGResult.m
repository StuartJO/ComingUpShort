function [FigureOutputLocs,LegendOutputLoc] = PlotFLaGResult(Features2Plot,PlotLabel,SAVEDIR,LgndSize)

load('GNM_FLaG_BestResults.mat','MdlBestFitIndv','EmpFit','Mdl_names','MdlBestFitAll')

FitName = {'max(\itKS\rm)','max(\itRMSE\rm)','max(\itr_d\rm )','\itTND','\itTF_{diff }','Degree correlation','Connection overlap (Jaccard)'};

cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
cmap2 = [0 0 0; cmap];

for i = 1:size(Features2Plot,1)
figure
FitStatIND = Features2Plot(i,1);
F1 = Features2Plot(i,2);
F2 = Features2Plot(i,3);

EmpDataF1 = EmpFit(:,F1);
EmpDataF2 = EmpFit(:,F2);

% MdlDataF1 = squeeze(MdlBestFitIndv(:,FitStatIND,F1,:));
% MdlDataF2 = squeeze(MdlBestFitIndv(:,FitStatIND,F2,:));
% MdlDataInd = repmat([2:11]',[1 size(MdlDataF2,2)]);

MdlDataF1 = [];
MdlDataF2 = [];
MdlDataInd = [];
for k = 1:10
    data2add = MdlBestFitAll{k,FitStatIND}(:,F1);
    MdlDataF1 = [MdlDataF1; data2add];
    MdlDataInd = [MdlDataInd; ones(size(data2add)).*(k+1)];
    data2add = MdlBestFitAll{k,FitStatIND}(:,F2);
    MdlDataF2 = [MdlDataF2; data2add];
end
% The order is randomised so the first points, which all belong to one
% model, aren't fully obscured by subsequent model data points
[~,randOrd] = sort(rand(length(MdlDataInd),1));
MdlDataF1 = MdlDataF1(randOrd);
MdlDataF2 = MdlDataF2(randOrd);
MdlDataInd = MdlDataInd(randOrd);

EmpInd = ones(size(EmpDataF1));

Xdata = [EmpDataF1(:); MdlDataF1(:)];
Ydata = [EmpDataF2(:); MdlDataF2(:)];

dataInd = [EmpInd(:); MdlDataInd(:)];

scatter(Xdata,Ydata,30,cmap2(dataInd,:),'filled','MarkerFaceAlpha',.1); %colormap(cmap2); clim([.5 11.5])

hold on
for k = 1:11
    meanValx = mean(Xdata(dataInd==k));
    meanValy = mean(Ydata(dataInd==k));

    if k == 1
        scatter(meanValx,meanValy,50,'filled','MarkerFaceColor','w','MarkerEdgeColor','k')
    
    else
        scatter(meanValx,meanValy,50,'filled','MarkerFaceColor',cmap2(k,:),'MarkerEdgeColor','k')
    
    end
end

if ismember(FitStatIND,[6 7])
title(['Networks with the best ',FitName{FitStatIND}])
else
title(['Networks with the lowest ',FitName{FitStatIND}])
end

xlabel(FitName{F1})
ylabel(FitName{F2})

set(gca,'FontSize',16)

annot = annotation(gcf, 'textbox',...
        [0,  .88, 0.0 0.0],...
        'String',PlotLabel{i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";

FigureOutputLocs{i} = [SAVEDIR,'/Fit',num2str(FitStatIND),'_Feat',num2str(F1),'_Feat',num2str(F2),'_Panel',PlotLabel{i},'.png'];
print(FigureOutputLocs{i},'-dpng','-r300')


end

Spos = gcf().Position;
TARGET_DPI = 300;
ScatterPlotWidth = round((TARGET_DPI/96)*Spos(3));
DesiredWidth = ScatterPlotWidth*LgndSize;

figure('Position',[50 100 1161 84])
clear ss
for i = 1:11
    hold on
    %ss(i) = scatter(-1,-1,100,cmap(i,:),'filled');
    ss(i) = plot(nan, nan,'Color','none','MarkerSize', 15,'Marker','o','MarkerFaceColor',cmap2(i,:));
end

leg = legend(ss,[{'Empirical'} Mdl_names],'Orientation','Horizontal','Location','northoutside','FontSize',16,'NumColumns',4,'Box','off');
leg.Position=[-0.0053    0.0349    1.0024    1.0024];
box off
axis off
%leg.Position=[-0.0056    0.2741    1.0001    0.4906];


Spos1 = gcf().Position;

SaveRes = round(96*(DesiredWidth/Spos1(3)));

LegendOutputLoc = [SAVEDIR,'/LEGEND_sz',num2str(LgndSize),'.png'];
print(LegendOutputLoc,'-dpng',['-r',num2str(SaveRes)])

