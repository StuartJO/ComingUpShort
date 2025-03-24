load('Hansen_networks.mat')

athr = triu2vec(adj{1},1);

A = adj{1};
Clu = clustering_coef_bu(A);

N = 200;
BCnormFactor = ((N-1)*(N-2))/2;
Bet = betweenness_bin(adj{1})';

Adeg = sum(adj{1},2);
NodeDist = sum(A_dist.*adj{1},2)./Adeg;

meanDeg = mean(Adeg);
meanClu = mean(Clu);
meanBet = mean(Bet);
meanNodeDist = mean(NodeDist);

Clo = closeness_bin(adj{1},1)';

A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = sum(A_dist.*A,2)./sum(A,2);
A_vals{5} = closeness_bin(A,1)';
A_vals{6} = mean(matching(A),2);

TF_A = corr([A_vals{1} A_vals{2} A_vals{3} A_vals{4} A_vals{5} A_vals{6}]);
GlobalA = [efficiency_bin(A) diffusion_efficiency(A) modularity_Q(A) transitivity_bu(A)];

MdlOptimFitType = zeros(100,6,10);
MdlFitType = zeros(100,6,10);
MdlDegCorr=zeros(100,10);
mdls=0:9;
NMdls=length(mdls);
cmap = [0.5 0.5 0.5; lines(7); 0.6941 0.3490 0.1569; [251,154,153]./255];
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
MdlName = ['Hansen_Timing_',num2str(TIMING),'_RandMdl',num2str(1),'_',AddMult,'_',LAW];
elseif mdl == 9
    MdlName = ['Hansen_Timing_',num2str(TIMING),'_TopoMdl_',num2str(1),'_',AddMult,'_',LAW];
else
MdlName = ['Hansen_Timing_',num2str(TIMING),'_Mdl',num2str(mdl),'_',AddMult,'_',LAW];
end
load([MdlName,'.mat'],'P','maxKS','b','optim_b','optim_maxKS','optim_DegCorr')
%MdlFits{mdlIND} = load([MdlName,'_TopoTypes.mat']);
load([MdlName,'_optim_FitType.mat'],'optim_FitType');
MdlFitType(:,:,mdlIND) = optim_FitType;
MdlDegCorr(:,mdlIND) = optim_DegCorr;
% MdlFits{mdlIND}.TND = pdist2(GlobalA,MdlFits{mdlIND}.TopoGlobal);
% MdlFits{mdlIND}.maxKS = maxKS;
% 
% MdlFitType(:,:,mdlIND) = [MdlFits{mdlIND}.maxKS max(MdlFits{mdlIND}.RMSE(:,[1 2 3 4]),[],2) 1-max(MdlFits{mdlIND}.TopogCorr,[],2) MdlFits{mdlIND}.TND' MdlFits{mdlIND}.TFdiff 1-MdlFits{mdlIND}.EdgeJaccard];
% 
% for i = 1:100
% B = zeros(200);
% B(optim_b{i})=1;
% B = B + B';
% [EdgeOverlap,GlobalTopo,RMSE,TopogCorr,TF_B] = TopoEvalMetrics(B,A,A_dist,A_vals);
% tfdiff = norm(TF_A-TF_B);
% tnd = pdist2(GlobalA,GlobalTopo);
% MdlOptimFitType(i,:,mdlIND) = [optim_maxKS(i) max(RMSE) 1-max(TopogCorr) tnd tfdiff 1-EdgeOverlap(:,4)];
% end

end

        end
    end
end

FtypeName = {'\itmax(KS)','\itmax(RMSE)','\itmax(r_d )','\itTND','\itTF_{diff }','Connection overlap (Jaccard)'};
% EmpFitType = {MdlFits.maxKS,max(MdlFits.RMSE(:,[1 2 3 4]),[],2),1-max(MdlFits.TopogCorr,[],2),MdlFits.TND,MdlFits.TF,1-MdlFits.EdgeJaccard};

EmpFits = load('C:\Users\Stuart\Documents\GitHub\PhysGenNetMdl\outputs\Schaefer_7net_iFOD2_acpc_lh_strThr_TopoComp.mat');

sub2use = [1:298 300:973];

EmpFitType(:,1) = triu2vec(EmpFits.maxKS(sub2use,sub2use,4),1);
EmpFitType(:,2) = triu2vec(squeeze(max(EmpFits.RMSE(sub2use,sub2use,4,1:4),[],4)),1);
EmpFitType(:,3) = 1-triu2vec(squeeze(max(EmpFits.TopogCorr(sub2use,sub2use,4,1:4),[],4)),1);
EmpFitType(:,4) = triu2vec(EmpFits.TopoDist(sub2use,sub2use,4),1);
EmpFitType(:,5) = triu2vec(EmpFits.TF(sub2use,sub2use,4),1);
EmpFitType(:,6) = triu2vec(EmpFits.EdgeOverlap{1}(sub2use,sub2use,4),1);

EmpDegCorr = triu2vec(EmpFits.DegCorr(sub2use,sub2use,4),1);

run = 0;
if run == 1

FIGLABEL = {'A','B','C','D','E'};

mkdir ./figures/emp

for i = 1:5
 cmap2 = [0 0 0; cmap];
MDLDATA = squeeze(MdlFitType(:,i,:));
MDLDATA_IND = ones(size(MDLDATA))+repmat(1:10,size(MDLDATA,1),1);
EMPDATA = EmpFitType(:,i);
EMPDATA_IND = ones(size(EMPDATA));

DATAi = [EMPDATA(:); MDLDATA(:)];
DATAind = [EMPDATA_IND(:); MDLDATA_IND(:)]; 
MDLDATA = 1-squeeze(MdlFitType(:,6,:));
EMPDATA = EmpFitType(:,6);
DATAj = [EMPDATA(:); MDLDATA(:)];
[PlotMain,PlotTop,PlotSide] = scatterWithKSden(DATAi,DATAj,DATAind,'grouping',DATAind,'grouping_colors',cmap2,'colormap',cmap2,'PlotColorbar','off','Annot',FIGLABEL{i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel(FtypeName{i})
ylabel('Connection overlap (Jaccard)')
print(['./figures/emp/Panel',FIGLABEL{i},'.png'],'-dpng','-r300')
end
% set(gcf, 'CurrentAxes', PlotTop)

FIGLABEL = {'F','G','H','I','J'};

for i = 1:5
 cmap2 = [0 0 0; cmap];
MDLDATA = squeeze(MdlFitType(:,i,:));
MDLDATA_IND = ones(size(MDLDATA))+repmat(1:10,size(MDLDATA,1),1);
EMPDATA = EmpFitType(:,i);
EMPDATA_IND = ones(size(EMPDATA));

DATAi = [EMPDATA(:); MDLDATA(:)];
DATAind = [EMPDATA_IND(:); MDLDATA_IND(:)]; 
MDLDATA = MdlDegCorr;
EMPDATA = EmpDegCorr;
DATAj = [EMPDATA(:); MDLDATA(:)];
[PlotMain,PlotTop,PlotSide] = scatterWithKSden(DATAi,DATAj,DATAind,'grouping',DATAind,'grouping_colors',cmap2,'colormap',cmap2,'PlotColorbar','off','Annot',FIGLABEL{i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel(FtypeName{i})
ylabel('Degree Correlation')
print(['./figures/emp/Panel',FIGLABEL{i},'.png'],'-dpng','-r300')
end


FIGLABEL = {'A','B','C','D','E'};

mkdir ./figures/emp

for i = 1:4
 cmap2 = [0 0 0; cmap];
MDLDATA = squeeze(MdlFitType(:,i+1,:));
MDLDATA_IND = ones(size(MDLDATA))+repmat(1:10,size(MDLDATA,1),1);
EMPDATA = EmpFitType(:,i+1);
EMPDATA_IND = ones(size(EMPDATA));

DATAi = [EMPDATA(:); MDLDATA(:)];
DATAind = [EMPDATA_IND(:); MDLDATA_IND(:)]; 
MDLDATA = squeeze(MdlFitType(:,1,:));
EMPDATA = EmpFitType(:,1);
DATAj = [EMPDATA(:); MDLDATA(:)];
[PlotMain,PlotTop,PlotSide] = scatterWithKSden(DATAi,DATAj,DATAind,'grouping',DATAind,'grouping_colors',cmap2,'colormap',cmap2,'PlotColorbar','off','Annot',FIGLABEL{i});
set(gcf, 'CurrentAxes', PlotMain)
xlabel(FtypeName{i+1})
ylabel('\itmax(KS)')
print(['./figures/emp/PanelS',FIGLABEL{i},'.png'],'-dpng','-r300')
end
end