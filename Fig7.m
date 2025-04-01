
% EdgeOverlap = zeros(100,10);
% 
% MdlFitType = zeros(100,6,10);
% MdlDegCorr=zeros(100,10);
mdls=1:10;
NMdls=length(mdls);
cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
 cmap2 = [0 0 0; cmap];
loadin = 1;
if loadin

BestFit = cell(10,5);
for i = 1:NMdls
MdlData = load(['EmpFLaG_Mdl',num2str(i),'_Add_exponential.mat'],'maxKS','maxRMSE','rd','TND','TFdiff','TFdiffnc','DegCorr','EdgeOverlap');
FIT{1} = MdlData.maxKS;
FIT{2} = MdlData.maxRMSE;
FIT{3} = MdlData.rd;
FIT{4} = MdlData.TND;
FIT{5} = MdlData.TFdiff;
EdgeJaccard = squeeze(MdlData.EdgeOverlap(:,:,4));
Nsub = size(EdgeJaccard,2);
for j = 1:5
    [~,I] = min(FIT{j});
    Iunique = unique(I);
    BestFit{i,j} = zeros(length(Iunique)*Nsub,7);
    for k = 1:5
        fitdata = FIT{k}(Iunique,:);
        BestFit{i,j}(:,k) = fitdata(:);
    end
    data1 = EdgeJaccard(Iunique,:);
    BestFit{i,j}(:,6) = data1(:);
    data2 = MdlData.DegCorr(Iunique,:);
    BestFit{i,j}(:,7) = data2(:);
end
end

FtypeName = {'max(\itKS\rm)','max(\itRMSE\rm)','max(\itr_d\rm )','\itTND','\itTF_{diff }','Connection overlap (Jaccard)'};
% EmpFitType = {MdlFits.maxKS,max(MdlFits.RMSE(:,[1 2 3 4]),[],2),1-max(MdlFits.TopogCorr,[],2),MdlFits.TND,MdlFits.TF,1-MdlFits.EdgeJaccard};

EmpFits = load('C:\Users\Stuart\Documents\GitHub\ComingUpShort\outputs\Scha400_7_lh_TopoMeasures_Thr70.mat');

sub2use = [1:298 300:973];

EmpFitType(:,1) = triu2vec(EmpFits.maxKS(sub2use,sub2use),1);
EmpFitType(:,2) = triu2vec(EmpFits.maxRMSE(sub2use,sub2use),1);
EmpFitType(:,3) = triu2vec(EmpFits.rd(sub2use,sub2use),1);
EmpFitType(:,4) = triu2vec(EmpFits.TND(sub2use,sub2use),1);
EmpFitType(:,5) = triu2vec(EmpFits.TFdiff(sub2use,sub2use),1);
EmpFitType(:,6) = triu2vec(EmpFits.EdgeOverlapJacc(sub2use,sub2use),1);

EmpDegCorr = triu2vec(EmpFits.DegCorr(sub2use,sub2use),1);

end

run = 1;
if run == 1

FIGLABEL = {'A','B','C','D','E'};

mkdir ./figures/emp

FIGLABEL = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
ITER = 1;

FIGLABEL_ITER = [1 1 2 3 4];

for i = 1:5
    for j = i
        figure
 cmap2 = [0 0 0; cmap];

    MDLDATA = [];
    MDLDATA_IND = [];
    for k = 1:10
        data2add = BestFit{k,i}(:,j);
        MDLDATA = [MDLDATA; data2add];
        MDLDATA_IND = [MDLDATA_IND; ones(size(data2add)).*(k+1)];
    end
     [~,randOrd] = sort(rand(length(MDLDATA_IND),1));
    MDLDATA = MDLDATA(randOrd);
    MDLDATA_IND = MDLDATA_IND(randOrd);


EMPDATA = EmpFitType(:,j);
EMPDATA_IND = ones(size(EMPDATA));

DATAi = [EMPDATA(:); MDLDATA(:)];
DATAind = [EMPDATA_IND(:); MDLDATA_IND(:)]; 

   MDLDATA = [];
    for k = 1:10
        data2add = BestFit{k,i}(:,6);
        MDLDATA = [MDLDATA; data2add];
   end
    MDLDATA = MDLDATA(randOrd);

EMPDATA = EmpFitType(:,6);
DATAj = [EMPDATA(:); MDLDATA(:)];
% [PlotMain,PlotTop,PlotSide] = scatterWithKSden(DATAi,DATAj,DATAind,'grouping',DATAind,'grouping_colors',cmap2,'colormap',cmap2,'PlotColorbar','off','Annot',FIGLABEL{i});
% set(gcf, 'CurrentAxes', PlotMain)
scatter(DATAi,DATAj,30,cmap2(DATAind,:),'filled','MarkerFaceAlpha',.1); colormap(cmap2); clim([.5 10.5])

hold on
for k = 1:11
    meanVali = mean(DATAi(DATAind==k));
    meanValj = mean(DATAj(DATAind==k));

if k == 1
scatter(meanVali,meanValj,50,'filled','MarkerFaceColor','w','MarkerEdgeColor','k')

else
scatter(meanVali,meanValj,50,'filled','MarkerFaceColor',cmap2(k,:),'MarkerEdgeColor','k')

end
end

xlabel(FtypeName{i})
ylabel('Connection overlap (Jaccard)')
if i ~= 1
title(['Networks with lowest ',FtypeName{i}])
end
set(gca,'FontSize',16)

annot = annotation(gcf, 'textbox',...
        [0,  .88, 0.0 0.0],...
        'String',FIGLABEL{FIGLABEL_ITER(i)},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
if i == 1
print(['./figures/emp/Fig7_',FIGLABEL{FIGLABEL_ITER(i)},'.png'],'-dpng','-r300')
else
print(['./figures/emp/PanelS_',FIGLABEL{FIGLABEL_ITER(i)},'.png'],'-dpng','-r300')
end
ITER = ITER+1;
close all
    end
end

FIGLABEL_ITER = [2 5 6 7 8];

for i = 1:5
    for j = i
        figure
 cmap2 = [0 0 0; cmap];

    MDLDATA = [];
    MDLDATA_IND = [];
    for k = 1:10
        data2add = BestFit{k,i}(:,j);
        MDLDATA = [MDLDATA; data2add];
        MDLDATA_IND = [MDLDATA_IND; ones(size(data2add)).*(k+1)];
    end
     [~,randOrd] = sort(rand(length(MDLDATA_IND),1));
    MDLDATA = MDLDATA(randOrd);
    MDLDATA_IND = MDLDATA_IND(randOrd);


EMPDATA = EmpFitType(:,j);
EMPDATA_IND = ones(size(EMPDATA));

DATAi = [EMPDATA(:); MDLDATA(:)];
DATAind = [EMPDATA_IND(:); MDLDATA_IND(:)]; 

   MDLDATA = [];
    for k = 1:10
        data2add = BestFit{k,i}(:,7);
        MDLDATA = [MDLDATA; data2add];
   end
    MDLDATA = MDLDATA(randOrd);

EMPDATA = EmpDegCorr;
DATAj = [EMPDATA(:); MDLDATA(:)];
% [PlotMain,PlotTop,PlotSide] = scatterWithKSden(DATAi,DATAj,DATAind,'grouping',DATAind,'grouping_colors',cmap2,'colormap',cmap2,'PlotColorbar','off','Annot',FIGLABEL{i});
% set(gcf, 'CurrentAxes', PlotMain)
scatter(DATAi,DATAj,30,cmap2(DATAind,:),'filled','MarkerFaceAlpha',.1); colormap(cmap2); clim([.5 10.5])
hold on
for k = 1:11
    meanVali = mean(DATAi(DATAind==k));
    meanValj = mean(DATAj(DATAind==k));

if k == 1
scatter(meanVali,meanValj,50,'filled','MarkerFaceColor','w','MarkerEdgeColor','k')

else
scatter(meanVali,meanValj,50,'filled','MarkerFaceColor',cmap2(k,:),'MarkerEdgeColor','k')

end
end

xlabel(FtypeName{i})
ylabel('Degree correlation')
%if i ~= 1
title(['Networks with lowest ',FtypeName{i}])
%end
set(gca,'FontSize',16)

annot = annotation(gcf, 'textbox',...
        [0,  .88, 0.0 0.0],...
        'String',FIGLABEL{FIGLABEL_ITER(i)},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";

if i == 1
print(['./figures/emp/Fig7_',FIGLABEL{FIGLABEL_ITER(i)},'.png'],'-dpng','-r300')
else
print(['./figures/emp/PanelS_',FIGLABEL{FIGLABEL_ITER(i)},'.png'],'-dpng','-r300')
end
ITER = ITER+1;
close all
    end
end


end

%%
runAnimation=0;

if runAnimation
EMPDATA = EmpFitType(:,1);
EMPDATA_IND = ones(size(EMPDATA));

EMPEDGEOVLP = EmpFitType(:,6);

scatter(EMPDATA,EMPEDGEOVLP,10,'k','filled','MarkerFaceAlpha',.1)

hold on

  MDLDATA = [];
    MDLDATA_IND = [];
    EdgeOverlap=[];
    MdlDegCorr = [];
    for k = 1:10
        data2add = BestFit{k,1}(:,1);
        MDLDATA = [MDLDATA; data2add];
        EdgeOverlap = [EdgeOverlap;BestFit{k,1}(:,6)];
        MdlDegCorr = [MdlDegCorr;BestFit{k,1}(:,7)];
        MDLDATA_IND = [MDLDATA_IND; ones(size(data2add)).*(k+1)];
    end
     [~,randOrd] = sort(rand(length(MDLDATA_IND),1));
    MDLDATA = MDLDATA(randOrd);
    EdgeOverlap = EdgeOverlap(randOrd);
    MdlDegCorr = MdlDegCorr(randOrd);
    MDLDATA_IND = MDLDATA_IND(randOrd);
for k = 1:10
meanMDLData(k) = mean(MDLDATA(MDLDATA_IND==(k+1)));
meanEdgeOverlap(k) = mean(EdgeOverlap(MDLDATA_IND==(k+1)));
meanMdlDegCorr(k) = mean(MdlDegCorr(MDLDATA_IND==(k+1)));
end
%MDLDATA = squeeze(MdlFitType(:,1,:));
%MDLDATA_IND = ones(size(MDLDATA))+repmat(1:10,size(MDLDATA,1),1);

scatter(MDLDATA(:),EdgeOverlap(:),10,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.5)
colormap(cmap2)



close all


figure('Position',[100  100 1920*.5 1080*.5])
Spos1 = gcf().Position;
ScatterPlotWidth = round((300/96)*1080);
DesiredWidth = ScatterPlotWidth;
SaveRes = round(96*(DesiredWidth/Spos1(3)));

scatter(EMPDATA,EMPEDGEOVLP,20,'k','filled','MarkerFaceAlpha',.1)
hold on
scatter(mean(EMPDATA),mean(EMPEDGEOVLP),50,'w','filled','MarkerEdgeColor','k');
start_xlim = xlim;
start_ylim = ylim;
scatter(MDLDATA(:),EdgeOverlap(:),20,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.15)
colormap(cmap2)
scatter(meanMDLData,meanEdgeOverlap,50,2:11,'filled','MarkerEdgeColor','k')

final_xlim = xlim;
final_ylim = ylim;

xlabel(FtypeName{1})
ylabel('Connection overlap')
set(gca,'FontSize',16)

R = linspace(0,1,60);
clim([1 11])
for t = 1:60
xlim(find_point_on_line(start_xlim,final_xlim,R(t)))
ylim(find_point_on_line(start_ylim,final_ylim,R(t)))

print(['./Animations/EmpKS/A0_',num2str(t),'.png'],'-dpng',['-r',num2str(SaveRes)])
end

figure('Position',[100  100 1920*.5 1080*.5])
scatter(EMPDATA,EmpDegCorr,20,'k','filled','MarkerFaceAlpha',.1)
hold on

scatter(mean(EMPDATA),mean(EmpDegCorr),50,'w','filled','MarkerEdgeColor','k');

start_xlim = xlim;
start_ylim = ylim;

scatter(MDLDATA(:),MdlDegCorr(:),20,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.15)
colormap(cmap2)

scatter(meanMDLData,meanMdlDegCorr,50,2:11,'filled','MarkerEdgeColor','k')

final_xlim = xlim;
final_ylim = ylim;

xlabel(FtypeName{1})
ylabel('Degree correlation')
set(gca,'FontSize',16)
clim([1 11])
print(['./Animations/EmpKS/DegCorr.png'],'-dpng',['-r',num2str(SaveRes)])

%%

L2 = [linspace(0,1,15) linspace(1,0,15)];

L = linspace(0,1,30);
for i = 1:4
figure('Position',[100  100 1920*.5 1080*.5])
EMPDATA = EmpFitType(:,i+1);

  MDLDATA = [];
    MDLDATA_IND = [];
    EdgeOverlap=[];
    MdlDegCorr = [];
    for k = 1:10
        data2add = BestFit{k,i+1}(:,i+1);
        MDLDATA = [MDLDATA; data2add];
        EdgeOverlap = [EdgeOverlap;BestFit{k,i+1}(:,6)];
        MdlDegCorr = [MdlDegCorr;BestFit{k,i+1}(:,7)];
        MDLDATA_IND = [MDLDATA_IND; ones(size(data2add)).*(k+1)];
    end
     [~,randOrd] = sort(rand(length(MDLDATA_IND),1));
    MDLDATA = MDLDATA(randOrd);
    EdgeOverlap = EdgeOverlap(randOrd);
    MdlDegCorr = MdlDegCorr(randOrd);
    MDLDATA_IND = MDLDATA_IND(randOrd);
for k = 1:10
meanMDLData(k) = mean(MDLDATA(MDLDATA_IND==(k+1)));
meanEdgeOverlap(k) = mean(EdgeOverlap(MDLDATA_IND==(k+1)));
meanMdlDegCorr(k) = mean(MdlDegCorr(MDLDATA_IND==(k+1)));
end


scatter(EMPDATA,EMPEDGEOVLP,20,'k','filled','MarkerFaceAlpha',.1);
hold on
scatter(mean(EMPDATA),mean(EMPEDGEOVLP),50,'w','filled','MarkerEdgeColor','k');
scatter(MDLDATA(:),EdgeOverlap(:),20,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.15)
colormap(cmap2)
clim([1 11])
scatter(meanMDLData,meanEdgeOverlap,50,2:11,'filled','MarkerEdgeColor','k')

xlabel(FtypeName{1+1})
ylabel('Connection overlap')
set(gca,'FontSize',16)

start_xlim = xlim;
start_ylim = ylim;

Xdata1_emp = EMPEDGEOVLP;
Xdata1_mdl = EdgeOverlap;

Xdata2_emp = EmpDegCorr;
Xdata2_mdl = MdlDegCorr;

for k = 1:10
meanXdata1_mdl(k) = mean(Xdata1_mdl(MDLDATA_IND==(k+1)));
meanXdata2_mdl(k) = mean(Xdata2_mdl(MDLDATA_IND==(k+1)));
end

clf

scatter(EMPDATA,Xdata2_emp,20,'k','filled','MarkerFaceAlpha',.1);
hold on
scatter(mean(EMPDATA),mean(Xdata2_emp),50,'w','filled','MarkerEdgeColor','k');
scatter(MDLDATA(:),Xdata2_mdl(:),20,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.15)
colormap(cmap2)
scatter(meanMDLData,meanXdata2_mdl,50,2:11,'filled','MarkerEdgeColor','k')
clim([1 11])
xlabel(FtypeName{i+1})
ylabel('Degree correlation')
set(gca,'FontSize',16)

final_xlim = xlim;
final_ylim = ylim;

for t = 1:30
clf

scatter(EMPDATA,find_point_on_line(Xdata1_emp,Xdata2_emp,L(t)),20,'k','filled','MarkerFaceAlpha',.1);
hold on

scatter(nanmean(EMPDATA),find_point_on_line(nanmean(Xdata1_emp),nanmean(Xdata2_emp),L(t)),50,'w','filled','MarkerEdgeColor','k');

scatter(MDLDATA(:),find_point_on_line(Xdata1_mdl(:),Xdata2_mdl(:),L(t)),20,MDLDATA_IND(:),'filled','MarkerFaceAlpha',.15)

colormap(cmap2)
scatter(meanMDLData,find_point_on_line(meanXdata1_mdl,meanXdata2_mdl,L(t)),50,2:11,'filled','MarkerEdgeColor','k')
clim([1 11])
xlabel(FtypeName{1+1})

if t < 15
ylabel('Connection overlap','Color',find_point_on_line([0 0 0],[1 1 1],L2(t)))
elseif t >=17
ylabel('Degree correlation','Color',find_point_on_line([0 0 0],[1 1 1],L2(t)))
else
ylabel('');
end


xlim(find_point_on_line(start_xlim,final_xlim,L(t)))
ylim(find_point_on_line(start_ylim,final_ylim,L(t)))

xlabel(FtypeName{i+1})
%ylabel('Degree correlation')
set(gca,'FontSize',16)

print(['./Animations/EmpKS/A',num2str(i),'_',num2str(t),'.png'],'-dpng',['-r',num2str(SaveRes)])


end


end

end