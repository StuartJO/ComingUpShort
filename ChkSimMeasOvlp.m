mdldata = load('Hansen_networks.mat');

A_dist = mdldata.A_dist;

adjs = mdldata.adj;
A = adjs{1};
d = triu2vec(mdldata.A_dist,1);

avec = triu2vec(adjs{1},1);

mdls = 0:9;

NMdls = length(mdls);

EvalMeasuresAll = cell(1,7);
for i = 1:7
    EvalMeasuresAll{i} = zeros(NMdls,10000);
end

MakePlot = 0;

A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = sum(A_dist.*A,2)./sum(A,2);
A_vals{5} = closeness_bin(A,1)';
A_vals{6} = mean(matching(A),2);

TF_A = corr([A_vals{1} A_vals{2} A_vals{3} A_vals{4} A_vals{5} A_vals{6}]);
GlobalA = [efficiency_bin(A) diffusion_efficiency(A) modularity_Q(A) transitivity_bu(A)];

for TIMING = 0
    for form = 1
        for law = 1

            if form == 0
                AddMult = 'Mult';
            else
            AddMult = 'Add';
            end
            if law == 0
                LAW = 'powerlaw';
            else
            LAW = 'exponential';
            end
            form_name = ['PhysMdls_Timing_',num2str(TIMING),'_',AddMult,'_',LAW];
for mdlIND = 1:NMdls
    mdl = mdls(mdlIND);
if mdl == 8
    LOADNAME = ['Hansen_Timing_',num2str(TIMING),'_RandMdl',num2str(1),'_',AddMult,'_',LAW];

elseif mdl == 9
    LOADNAME = ['Hansen_Timing_',num2str(TIMING),'_TopoMdl_',num2str(1),'_',AddMult,'_',LAW];
else
LOADNAME = ['Hansen_Timing_',num2str(TIMING),'_Mdl',num2str(mdl),'_',AddMult,'_',LAW];
end
Output = load([LOADNAME,'.mat']);
FitTypes = load([LOADNAME,'_TopoTypes.mat']);

FtypeName = {'\itmax(KS)','\itmax(RMSE)','\itmax(r_d )','\itTND','\itTF_{diff }','Connection recovery','Degree correlation'};

    EvalMeasuresAll{1}(mdlIND,:) = Output.maxKS;
    EvalMeasuresAll{2}(mdlIND,:) = mean(FitTypes.RMSE(:,1:4),2);
    EvalMeasuresAll{3}(mdlIND,:) = 1-mean(FitTypes.TopogCorr(:,1:4),2);
    EvalMeasuresAll{4}(mdlIND,:) = pdist2(FitTypes.TopoGlobal,GlobalA);
    EvalMeasuresAll{5}(mdlIND,:) = FitTypes.TFdiff;

    EvalMeasuresAll{6}(mdlIND,:) = FitTypes.EdgeOverlap(:,1);
    EvalMeasuresAll{7}(mdlIND,:) = Output.DegCorr(:,1);


MAP_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

if MakePlot

figure('Position',[1  49 1382 747])
subPos = zeros(25,4);
for i = 1:5
    for j = 1:5
        subplot(5,5,sub2ind([5 5],i,j))
        if i == j 
        histogram(EvalMeasuresAll{i}(mdlIND,:))
        xlabel(FtypeName{i})
        elseif i > j
        scatter(EvalMeasuresAll{i}(mdlIND,:),EvalMeasuresAll{j}(mdlIND,:),30,EvalMeasuresAll{6}(mdlIND,:),'filled')
        colormap(gca,"parula")
        xlabel(FtypeName{i})
        ylabel(FtypeName{j})
        elseif i < j
        scatter(EvalMeasuresAll{i}(mdlIND,:),EvalMeasuresAll{j}(mdlIND,:),30,EvalMeasuresAll{7}(mdlIND,:),'filled')
        colormap(gca,"sky")
        xlabel(FtypeName{i})
        ylabel(FtypeName{j})
        end
        subPos(sub2ind([5 5],i,j),:) = get(gca,'Position');
        if i==4 && j == 5
            c1 = colorbar;            
            c1.Label.String = 'Degree correlation';
            c1.FontSize = 16;
        elseif i==5 && j == 4
            c2 = colorbar;
            
            c2.Label.String = 'Connection recovery';
            c2.FontSize = 16;
        end
    end
end
c1.Position = [.92 subPos(21,2) 0.015 .38];
c2.Position = [.92 (subPos(1,2)+subPos(1,4))-(.38) 0.015 .38];
exportgraphics(gcf,['CompFits',MAP_names{mdlIND},'.png'],'Resolution',300)
end

end
        end
    end
end

ANNOTS = {'A','B','C','D','E','F','G','H','I','J'};

Thr = [1 2.5 5 10];
Ovlp = zeros(5,5,4,10);
for i = 1:5
    for j = 1:5
        for k = 1:10
x = EvalMeasuresAll{i}(k,:);
y = EvalMeasuresAll{j}(k,:);
ISNAN = isnan(x)|isnan(y);
x(ISNAN)=[];
y(ISNAN)=[];
MeasureCorr{k}(i,j) = corr(x',y','type','Spearman');
[xsort,xsortind] = sort(x,'ascend');
[ysort,ysortind] = sort(y,'ascend');
for TopThrInd = 1:length(Thr)
    cutoff = floor(Thr(TopThrInd)/100 * length(x));
    Ovlp(i,j,TopThrInd,k) = length(intersect(xsortind(1:cutoff),ysortind(1:cutoff)))/cutoff;
end

        end
    end
end

FtypeName = {'\itmax(KS)','\itmax(RMSE)','\itmax(r_d )','\itTND','\itTF_{diff }','Connection recovery','Degree correlation'};

cmap_orig = cmocean('Balance',300);
cmap = cmap_orig(33:268,:);
%tiledlayout(2,5)
for i = 1:10
    %subplot(2,5,i)
    %nexttile
    figure
    imagesc(MeasureCorr{i})
    xticks(1:5)
    yticks(1:5)
    xticklabels(FtypeName(1:5))
    yticklabels(FtypeName(1:5))
    xtickangle(45)
    axis square
    clim([-1 1])
    colormap(cmap)
    title(MAP_names{i})
    c = colorbar;
    c.Label.String = 'Spearman correlation';
    c.FontSize = 16;
    c.LineWidth = 1.5;
    set(gca, 'LineWidth',1.5)
    set(gca,'FontSize',16)

    annot = annotation(gcf, 'textbox',...
        [0,  .9, 0.0 0.0],...
        'String',ANNOTS{i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
    set(gca,'TickLength',[.025 .1])
    %set(gca, 'TickDir', 'out') 
    print(['FitTypeCorr_',ANNOTS{i},'.png'],'-dpng')
end

for k = 1:4
for i = 1:10
    %subplot(2,5,i)
    %nexttile
    figure
    imagesc(squeeze(Ovlp(:,:,k,i)))
    xticks(1:5)
    yticks(1:5)
    xticklabels(FtypeName(1:5))
    yticklabels(FtypeName(1:5))
    xtickangle(45)
    axis square
    clim([0 1])
    colormap(parula)
    title(MAP_names{i})
    c = colorbar;
    c.Label.String = ['Prop. same top ',num2str(Thr(k)),'% networks'];
    c.FontSize = 16;
    c.LineWidth = 1.5;
    set(gca, 'LineWidth',1.5)
    set(gca,'FontSize',16)

    annot = annotation(gcf, 'textbox',...
        [0,  .9, 0.0 0.0],...
        'String',ANNOTS{i},...
        'LineStyle','none',...
        'FitBoxToText','on',...
        'FontSize', 28, ...
        'FitBoxToText','off');
    annot.VerticalAlignment = "bottom";
    set(gca,'TickLength',[.025 .1])
    %set(gca, 'TickDir', 'out') 
    print(['FitTypeTop',num2str(k),'_',ANNOTS{i},'.png'],'-dpng')
end
end

for i = 1:10
CorrAll(:,:,i) = MeasureCorr{i};
CorrAlltriu(i,:) = triu2vec(MeasureCorr{i});
end

for k = 1:5
    data = [];
    for i = 1:10
        data(i,:) = triu2vec(squeeze(CorrAll(k,:,i)));
    end
MeanCorrAll(k) = mean(data(:));
end
% 
% 
% iter = 1;
% for i = 1:10
%     for j = 1:4
%         data = squeeze(Ovlp(:,:,j,i));
%         subplot(4,10,iter)
%         imagesc(data)
%         clim([0 1])
%         iter = iter+1;
%         axis square
%         meanOvlp(j,i) = mean(triu2vec(data,1));
%     end
% end
% 
% iter = 1;
% for i = 1:5
%     for j = 1:5
%         for k = 1
% x = EvalMeasuresAll{i}(k,:);
% y = EvalMeasuresAll{j}(k,:);
% subplot(5,5,iter)
% [xsort,xsortind] = sort(x,'ascend');
% [ysort,ysortind] = sort(y,'ascend');
% 
% Intop = zeros(length(ysortind),1);
% cutoff = Thr(1)/100 * length(x);
% 
% Intop(xsortind(1:cutoff))=1;
% Intop(ysortind(1:cutoff))=Intop(ysortind(1:cutoff))+2;
% 
% scatter(x,y,30,Intop,'filled')
% clim([0 3])
% iter = iter+1;
%         end
%     end
% end
% 
% cmap = cmocean('Balance',256);