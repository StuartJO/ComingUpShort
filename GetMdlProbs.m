mdldata = load('Hansen_networks.mat');

mkdir .\figures\ProbPlot

A_dist = mdldata.A_dist;

A = mdldata.adj{1};

d = triu2vec(mdldata.A_dist,1);

avec = triu2vec(A,1);

mdls = 0:9;

NMdls = length(mdls);

prob=zeros(NMdls,19900);
ITER = 1;
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
            form_name = ['PhysMdls_',AddMult,'_',LAW];
for mdlIND = 1:NMdls
    mdl = mdls(mdlIND);
if mdl == 9
    Output = load(['GNM_TopoMdl',num2str(2),'_',AddMult,'_',LAW,'.mat']);
else
Output = load(['GNM_Mdl',num2str(mdl),'_',AddMult,'_',LAW,'.mat']);
end

prob(mdlIND,:) = Output.BestFit.maxKS.meanProb;

end

%mkdir(['.\figures\ProbPlot\',form_name])


MAP_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

Mdls2Use = [8 4 9 2];

if form == 1 && law == 1
ANNOTS = {'A','B','C','D'};
else
    ANNOTS_ = [{'A','B','C','D'};{'E','F','G','H'};{'I','J','K','L'}];
    ANNOTS = ANNOTS_(ITER,:);
ITER = ITER+1;
end

for i = 1:4
mdl = Mdls2Use(i);
dvec= triu2vec(A_dist);
probs = prob(mdl,:);
avec = triu2vec(A);
[~,sortOrd] = sort(avec,'ascend');

[PlotMain,PlotTop,PlotSide] = scatterWithKSden(dvec(sortOrd),probs(sortOrd),avec(sortOrd),'Yscale','log','grouping',avec(sortOrd),'grouping_colors',[0 0 0; 1 0 0],'colormap',[0 0 0; 1 0 0],'PlotColorbar','off','Annot',ANNOTS{i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel('Distance (mm)')
ylabel('Probability ({\itP_{ij}})','Position',[-9.756064875641492,0,-1])
set(gcf, 'CurrentAxes', PlotTop)
title(MAP_names{mdl})
if form == 1 && law == 1
    print(['./figures/ProbPlot/Main',ANNOTS{i},'.png'],'-dpng','-r300')
else
    print(['./figures/ProbPlot/Supp',ANNOTS{i},'.png'],'-dpng','-r300')
end
end

if form == 1 && law == 1
ANNOTS = {'A','B','C','D','E','F','G','H','I','J'};


for i = 1:10
mdl = i;
dvec= triu2vec(A_dist);
probs = prob(mdl,:);
avec = triu2vec(A);
[~,sortOrd] = sort(avec,'ascend');

[PlotMain,PlotTop,PlotSide] = scatterWithKSden(dvec(sortOrd),probs(sortOrd),avec(sortOrd),'Yscale','log','grouping',avec(sortOrd),'grouping_colors',[0 0 0; 1 0 0],'colormap',[0 0 0; 1 0 0],'PlotColorbar','off','Annot',ANNOTS{i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel('Distance (mm)')
ylabel('Probability ({\itP_{ij}})')
set(gcf, 'CurrentAxes', PlotTop)
title(MAP_names{mdl})
print(['./figures/ProbPlot/SuppALL_',ANNOTS{i},'.png'],'-dpng','-r300')
end

end

        end
    end

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ( ./figures/ProbPlot/SuppA.png ./figures/ProbPlot/SuppB.png ./figures/ProbPlot/SuppC.png ./figures/ProbPlot/SuppD.png +append ) ( ./figures/ProbPlot/SuppE.png ./figures/ProbPlot/SuppF.png ./figures/ProbPlot/SuppG.png ./figures/ProbPlot/SuppH.png +append ) ( ./figures/ProbPlot/SuppI.png ./figures/ProbPlot/SuppJ.png ./figures/ProbPlot/SuppK.png ./figures/ProbPlot/SuppL.png +append ) -append ./figures/Supp.png']);

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ( ./figures/ProbPlot/SuppALL_A.png ./figures/ProbPlot/SuppALL_B.png ./figures/ProbPlot/SuppALL_C.png ./figures/ProbPlot/SuppALL_D.png ./figures/ProbPlot/SuppALL_E.png +append ) ( ./figures/ProbPlot/SuppALL_F.png ./figures/ProbPlot/SuppALL_G.png ./figures/ProbPlot/SuppALL_H.png ./figures/ProbPlot/SuppALL_I.png ./figures/ProbPlot/SuppALL_J.png +append ) -append ./figures/SuppALL.png']);

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ./figures/ProbPlot/MainA.png ./figures/ProbPlot/MainB.png ./figures/ProbPlot/MainC.png ./figures/ProbPlot/MainD.png +append ./figures/Main.png']);
