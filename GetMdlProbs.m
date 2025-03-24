mdldata = load('Hansen_networks.mat');

mkdir .\figures\ProbPlot

A_dist = mdldata.A_dist;

A = mdldata.adj{1};

d = triu2vec(mdldata.A_dist,1);

avec = triu2vec(A,1);

dist_thr{2} = d<30;
dist_thr{3} = d>=30 & d<=90;
dist_thr{4}= d>90;
dist_thr{1}= d>=0;


d_short = d<30;
d_mid = d>=30 & d<=90;
d_long = d>90;

mdls = 0:9;

NMdls = length(mdls);

MdlFit = zeros(NMdls,100,1);
DegCorr = zeros(NMdls,100,1);
EdgeOverlap = zeros(NMdls,100,4);
EdgeProp = zeros(NMdls,100,4);
Edge_cdf = zeros(100,length(0:160),NMdls);

MdlEdges = zeros(NMdls,length(d));

f = zeros(NMdls,length(d));
Eta = zeros(NMdls,1);
Gam = zeros(NMdls,1);
Alpha = zeros(NMdls,1);

MdlDeg = zeros(NMdls,200);

prob=zeros(NMdls,19900);
ITER = 1;
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
Output = load(['Hansen_Timing_',num2str(TIMING),'_RandMdl',num2str(1),'_',AddMult,'_',LAW,'.mat']);
elseif mdl == 9
    Output = load(['Hansen_Timing_',num2str(TIMING),'_TopoMdl_',num2str(1),'_',AddMult,'_',LAW,'.mat']);
else
Output = load(['Hansen_Timing_',num2str(TIMING),'_Mdl',num2str(mdl),'_',AddMult,'_',LAW,'.mat']);
end

prob(mdlIND,:) = Output.optim_meanProb;

end


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
ylabel('Probability ({\itP_{ij}})')
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
end

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ( ./figures/ProbPlot/SuppA.png ./figures/ProbPlot/SuppB.png ./figures/ProbPlot/SuppC.png ./figures/ProbPlot/SuppD.png +append ) ( ./figures/ProbPlot/SuppE.png ./figures/ProbPlot/SuppF.png ./figures/ProbPlot/SuppG.png ./figures/ProbPlot/SuppH.png +append ) ( ./figures/ProbPlot/SuppI.png ./figures/ProbPlot/SuppJ.png ./figures/ProbPlot/SuppK.png ./figures/ProbPlot/SuppL.png +append ) -append ./figures/Supp.png']);

system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ( ./figures/ProbPlot/SuppALL_A.png ./figures/ProbPlot/SuppALL_B.png ./figures/ProbPlot/SuppALL_C.png ./figures/ProbPlot/SuppALL_D.png ./figures/ProbPlot/SuppALL_E.png +append ) ( ./figures/ProbPlot/SuppALL_F.png ./figures/ProbPlot/SuppALL_G.png ./figures/ProbPlot/SuppALL_H.png ./figures/ProbPlot/SuppALL_I.png ./figures/ProbPlot/SuppALL_J.png +append ) -append ./figures/SuppALL.png']);


system(['"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe" convert ./figures/ProbPlot/MainA.png ./figures/ProbPlot/MainB.png ./figures/ProbPlot/MainC.png ./figures/ProbPlot/MainD.png +append ./figures/Main.png']);
