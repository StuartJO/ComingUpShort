
ImgMagLoc = '"C:\Program Files\ImageMagick-7.1.1-Q16-HDRI/magick.exe"';

%%

[GNMmaxKS{1},LgdLoc] = PlotGNMPerformance('BestMdls_GNM_maxKS_Add_exponential.mat',{'A','B','C','D',''});

system([ImgMagLoc,' montage ',GNMmaxKS{1}{1},' ',GNMmaxKS{1}{2},' ',GNMmaxKS{1}{3},' -tile 3x1 -geometry +0+0 miff:- | magick convert - ',GNMmaxKS{1}{4},' ',LgdLoc,' -append ./figures/Figure2.png']);

system([ImgMagLoc,' convert ',GNMmaxKS{1}{5},' ',LgdLoc,' -append ./figures/FigureS5.png']);

close all

GNMmaxKS{2} = PlotGNMPerformance('BestMdls_GNM_maxKS_Add_powerlaw.mat',{'A','B','C','A','A'});
GNMmaxKS{3} = PlotGNMPerformance('BestMdls_GNM_maxKS_Mult_exponential.mat',{'D','E','F','B','B'});
GNMmaxKS{4} = PlotGNMPerformance('BestMdls_GNM_maxKS_Mult_powerlaw.mat',{'G','H','I','C','C'});
close all
FIGURES2ADD = [GNMmaxKS{2}{1},' ',GNMmaxKS{2}{2},' ',GNMmaxKS{2}{3},' ',GNMmaxKS{3}{1},' ',GNMmaxKS{3}{2},' ',GNMmaxKS{3}{3},' ',GNMmaxKS{4}{1},' ',GNMmaxKS{4}{2},' ',GNMmaxKS{4}{3}];

system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 3x3 -geometry +0+0 miff:- | magick convert - ',LgdLoc,' -append ./figures/FigureS1.png']);

FIGURES2ADD = [GNMmaxKS{2}{4},' ',GNMmaxKS{3}{4},' ',GNMmaxKS{4}{4},' ',LgdLoc];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 1x4 -geometry +0+0 ./figures/FigureS3.png']);

close all

[GNMTopomaxKS{1},LgdLocTopo] = PlotGNMPerformance('BestMdls_TopoMdl_GNM_maxKS_Mult_powerlaw.mat',{'A','B','C','A','A'});
GNMTopomaxKS{2} = PlotGNMPerformance('BestMdls_TopoMdl_GNM_maxKS_Mult_exponential.mat',{'D','E','F','B','B'});
GNMTopomaxKS{3} = PlotGNMPerformance('BestMdls_TopoMdl_GNM_maxKS_Add_powerlaw.mat',{'G','H','I','C','C'});
GNMTopomaxKS{4} = PlotGNMPerformance('BestMdls_TopoMdl_GNM_maxKS_Add_exponential.mat',{'J','K','L','D','D'});

FIGURES2ADD = [GNMTopomaxKS{1}{1},' ',GNMTopomaxKS{1}{2},' ',GNMTopomaxKS{1}{3},' ',GNMTopomaxKS{2}{1},' ',GNMTopomaxKS{2}{2},' ',GNMTopomaxKS{2}{3},' ',GNMTopomaxKS{3}{1},' ',GNMTopomaxKS{3}{2},' ',GNMTopomaxKS{3}{3},' ',GNMTopomaxKS{4}{1},' ',GNMTopomaxKS{4}{2},' ',GNMTopomaxKS{4}{3}];

system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 3x4 -geometry +0+0 miff:- | magick convert - ',LgdLocTopo,' -append ./figures/FigureS2.png']);

FIGURES2ADD = [GNMTopomaxKS{1}{4},' ',GNMTopomaxKS{2}{4},' ',GNMTopomaxKS{3}{4},' ',GNMTopomaxKS{4}{4},' ',LgdLocTopo];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 1x5 -geometry +0+0 ./figures/FigureS4.png']);

close all

GNMDegCorr{1} = PlotGNMPerformance('BestMdls_GNM_DegCorr_Add_exponential.mat',{'A','B','C','A','A'});
GNMDegCorr{2} = PlotGNMPerformance('BestMdls_GNM_DegCorr_Add_powerlaw.mat',{'D','E','F','B','B'});
GNMDegCorr{3} = PlotGNMPerformance('BestMdls_GNM_DegCorr_Mult_exponential.mat',{'G','H','I','C','C'});
GNMDegCorr{4} = PlotGNMPerformance('BestMdls_GNM_DegCorr_Mult_powerlaw.mat',{'J','K','L','D','D'});

FIGURES2ADD = [GNMDegCorr{1}{1},' ',GNMDegCorr{1}{2},' ',GNMDegCorr{1}{3},' ',GNMDegCorr{2}{1},' ',GNMDegCorr{2}{2},' ',GNMDegCorr{2}{3},' ',GNMDegCorr{3}{1},' ',GNMDegCorr{3}{2},' ',GNMDegCorr{3}{3},' ',GNMDegCorr{4}{1},' ',GNMDegCorr{4}{2},' ',GNMDegCorr{4}{3}];

system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 3x4 -geometry +0+0 miff:- | magick convert - ',LgdLocTopo,' -append ./figures/FigureS7.png']);

FIGURES2ADD = [GNMDegCorr{1}{4},' ',GNMDegCorr{2}{4},' ',GNMDegCorr{3}{4},' ',GNMDegCorr{4}{4},' ',LgdLoc];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 1x5 -geometry +0+0 ./figures/FigureS8.png']);

close all

GNM_WB_maxKS = PlotGNMPerformance('BestMdls_WB_GNM_maxKS_Add_exponential.mat',{'A','B','C','D','A','B','C','A','B','C'});

system([ImgMagLoc,' montage ',GNM_WB_maxKS{1},' ',GNM_WB_maxKS{2},' ',GNM_WB_maxKS{3},' -tile 3x1 -geometry +0+0 miff:- | magick convert - ',GNM_WB_maxKS{4},' ',LgdLoc,' -append ./figures/FigureS9.png']);

FIGURES2ADD = [GNM_WB_maxKS{5},' ',GNM_WB_maxKS{6},' ',GNM_WB_maxKS{7},' ',LgdLoc];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 1x4 -geometry +0+0 ./figures/FigureS10.png']);

close all


PlotLabels = {'A','B','C','D','E','F','G','A','H','I';'A','B','C','D','E','F','G','C','H','I';'A','B','C','D','E','F','G','B','H','I'};
FigureOutputLocs = PlotGNMProbabilities('BestMdls_GNM_maxKS_Add_exponential.mat',PlotLabels);
FIGURES2ADD = [FigureOutputLocs{1,8},' ',FigureOutputLocs{3,8}];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 2x1 -geometry +0+0 ./figures/Figure4.png']);

FIGURES2ADD = [FigureOutputLocs{1,1},' ',FigureOutputLocs{1,2},' ',FigureOutputLocs{1,3},' ',FigureOutputLocs{1,4},' ',FigureOutputLocs{1,5},' ',FigureOutputLocs{1,6},' ',FigureOutputLocs{1,7},' ',FigureOutputLocs{1,9},' ',FigureOutputLocs{1,10}];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 3x3 -geometry +0+0 ./figures/FigureS11.png']);

FIGURES2ADD = [FigureOutputLocs{3,1},' ',FigureOutputLocs{3,2},' ',FigureOutputLocs{3,3},' ',FigureOutputLocs{3,4},' ',FigureOutputLocs{3,5},' ',FigureOutputLocs{3,6},' ',FigureOutputLocs{3,7},' ',FigureOutputLocs{3,9},' ',FigureOutputLocs{3,10}];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 3x3 -geometry +0+0 ./figures/FigureS12.png']);

PlotLabels = {'A','B','C','D','E','F','G','H','I','J';'A','B','C','D','E','F','G','H','I','J';'A','B','C','D','E','F','G','H','I','J'};
FigureOutputLocs = PlotGNMProbabilities('BestMdls_WB_GNM_maxKS_Add_exponential.mat',PlotLabels);

%%

PlotIterativeRewirings

close all

%%

PlotEmpiricalKS('maxKS','edge','parc',100:100:1000,'iFOD2')
print('./figures/EmpParcAll_A.png','-dpng','-r300')
close all
PlotEmpiricalKS('maxKS','DegCorr','parc',100:100:1000,'iFOD2')
print('./figures/EmpParcAll_B.png','-dpng','-r300')
close all
PlotEmpiricalKS('edge','DegCorr','parc',100:100:1000,'iFOD2')
print('./figures/EmpParcAll_C.png','-dpng','-r300')
close all
system([ImgMagLoc,' montage ./figures/EmpParcAll_A.png ./figures/EmpParcAll_B.png ./figures/EmpParcAll_C.png -tile 1x3 -geometry +0+0 ./figures/Figure6.png']);

PlotEmpiricalKS('maxKS','edge','parc',100:100:1000,'FACT')
print('./figures/FACTEmpParcAll_A.png','-dpng','-r300')
close all
PlotEmpiricalKS('maxKS','DegCorr','parc',100:100:1000,'FACT')
print('./figures/FACTEmpParcAll_B.png','-dpng','-r300')
close all
PlotEmpiricalKS('edge','DegCorr','parc',100:100:1000,'FACT')
print('./figures/FACTEmpParcAll_C.png','-dpng','-r300')
close all
system([ImgMagLoc,' montage ./figures/FACTEmpParcAll_A.png ./figures/FACTEmpParcAll_B.png ./figures/FACTEmpParcAll_C.png -tile 1x3 -geometry +0+0 ./figures/FigureS14.png']);

%%
mkdir './figures/FLaGResult'
PlotLabels = {'A','B','C','D','E','F','G','H','I','J'};
Features2Plot = [1 1 6;1 1 7];
[FLAG1,FLAG1_leg] = PlotFLaGResult(Features2Plot,PlotLabels,'./figures/FLaGResult',2);

system([ImgMagLoc,' montage ',FLAG1{1},' ',FLAG1{2},' -tile 2x1 -geometry +0+0 miff:- | magick convert - ',FLAG1_leg,' -append ./figures/Figure7.png']);

Features2Plot = zeros(8,3);
for i = 2:5
Features2Plot(i-1,:) = [i i 6];
Features2Plot(i+3,:) = [i i 7];
end
[FLAG2,FLAG2_leg] = PlotFLaGResult(Features2Plot,PlotLabels,'./figures/FLaGResult',4);

FIGURES2ADD = [FLAG2{1},' ',FLAG2{2},' ',FLAG2{3},' ',FLAG2{4},' ',FLAG2{5},' ',FLAG2{6},' ',FLAG2{7},' ',FLAG2{8}];
system([ImgMagLoc,' montage ',FIGURES2ADD,' -tile 4x2 -geometry +0+0 miff:- | magick convert - ',FLAG2_leg,' -append ./figures/Figure8.png']);

%%