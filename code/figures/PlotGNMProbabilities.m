function FigureOutputLocs = PlotGNMProbabilities(INPUT,PlotLabels,SAVEDIR)

load(INPUT,'meanProb','Mdl_names')

NMdls = length(Mdl_names);

FigureOutputLocs = cell(size(PlotLabels));

if nargin < 3
    [~,name] = fileparts(INPUT);
    SAVEDIR=name;
end

if size(meanProb,2)==19900
    mdldata = load('Hansen_networks.mat');
else
    mdldata = load('Hansen_networks_WB.mat');
end
A_dist = mdldata.A_dist;
A = mdldata.adj{1};
dvec= triu2vec(A_dist);
avec = triu2vec(A,1);


for i = 1:NMdls
probs = meanProb(i,:);
[~,sortOrd] = sort(avec,'ascend');

[PlotMain,PlotTop,PlotSide] = scatterWithKSden(dvec(sortOrd),probs(sortOrd),avec(sortOrd),'Yscale','log','grouping',avec(sortOrd),'grouping_colors',[0 0 0; 1 0 0],'colormap',[0 0 0; 1 0 0],'PlotColorbar','off','Annot',PlotLabels{1,i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel('Distance (mm)')
ylabel('Probability ({\itP_{ij}})')
set(gcf, 'CurrentAxes', PlotTop)
title(Mdl_names{i})


FigureOutputLocs{1,i} = ['./figures/',SAVEDIR,'/ProbsScatter',num2str(i),'.png'];
print(FigureOutputLocs{1,i},'-dpng','-r300')
end

ptiles = [5 10 25 50 100];

maxd = round(max(dvec));
TopProb = zeros(maxd,length(ptiles),10);

EXISTprob = zeros(NMdls,maxd);
NOEXISITprob = zeros(NMdls,maxd);

for j = 1:maxd
dthr = j;
for i = 1:NMdls

    dthrlog = dvec>=dthr-5 & dvec<=dthr+5;
    athr = avec(dthrlog);

EXISTprob(i,j) = median(meanProb(i,dthrlog & avec==1));
NOEXISITprob(i,j) = median(meanProb(i,dthrlog & avec==0));

    probthr = meanProb(i,dthrlog);

    for k = 1:length(ptiles)
        ptilethr = prctile(probthr,100-ptiles(k));
        TopProb(j,k,i) = mean(athr(probthr>=ptilethr));
    end

end
end

ptile_leg = cell(1,length(ptiles));
for i = 1:length(ptiles)
ptile_leg{i} = [num2str(ptiles(i)),'%'];
end
ptile_leg{i} = 'Null';
%figure('Position',[1 1 1382 747])
for i = 1:NMdls
    %subplot(2,5,i)
    figure('Position',[100   100   706   694])
plot(EXISTprob(i,:),'k','LineWidth',3)
hold on
plot(NOEXISITprob(i,:),'r','LineWidth',3)
set(gca,'Yscale','log')
xlabel('Distance (mm)')
ylabel('Median probability ({\itP_{ij}})')
title(Mdl_names{i})
set(gca,'FontSize',16)
set(gca, 'LineWidth',1.5)
annot = annotation(gcf, 'textbox',...
        [0,  .92, 0.0 0.0],...
        'String',PlotLabels{2,i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 36, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";

FigureOutputLocs{2,i} = ['./figures/',SAVEDIR,'/MedianProbMdl_',num2str(i),'.png'];
print(FigureOutputLocs{2,i},'-dpng','-r300')
end

sky_cmap = [sky(length(ptiles)-1); 0 0 0];
figure('Position',[1 1 1382 747])
for i = 1:NMdls
    %subplot(2,5,i)
    figure('Position',[100   100   706   694])
    for j = 1:length(ptiles)
plot(squeeze(TopProb(:,j,i)),'Color',sky_cmap(j,:),'LineWidth',2)
hold on
    end
legend(ptile_leg,'FontSize',14)

xlabel('Distance (mm)')
ylabel({'Proportion of empirical','connections in top X%'})
title(Mdl_names{i})
set(gca,'FontSize',16)
set(gca, 'LineWidth',1.5)
annot = annotation(gcf, 'textbox',...
        [0,  .92, 0.0 0.0],...
        'String',PlotLabels{3,i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 36, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";

FigureOutputLocs{3,i} = ['./figures/',SAVEDIR,'/ProbTopMdl_',num2str(i),'.png'];
print(FigureOutputLocs{3,i},'-dpng','-r300')

end
close all