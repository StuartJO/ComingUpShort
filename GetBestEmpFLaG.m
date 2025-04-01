BestFit = cell(10,7);
for i = 1:10
MdlData = load(['./outputs/EmpFLaG_Mdl',num2str(i),'_Add_exponential.mat'],'maxKS','maxRMSE','rd','TND','TFdiff','TFdiffnc','DegCorr','EdgeOverlap');
FIT{1} = MdlData.maxKS;
FIT{2} = MdlData.maxRMSE;
FIT{3} = MdlData.rd;
FIT{4} = MdlData.TND;
FIT{5} = MdlData.TFdiff;
FIT{6} = MdlData.DegCorr;
FIT{7} = squeeze(MdlData.EdgeOverlap(:,:,4));
Nsub = size(FIT{7},2);
for j = 1:7
    [~,I] = min(FIT{j});
    Iunique = unique(I);
    BestFit{i,j} = zeros(length(Iunique)*Nsub,7);
    for k = 1:7
        fitdata = FIT{k}(Iunique,:);
        BestFit{i,j}(:,k) = fitdata(:);
    end
end
end

FtypeName = {'max(\itKS\rm)','max(\itRMSE\rm)','max(\itr_d\rm )','\itTND','\itTF_{diff }','Degree correlation','Connection overlap (Jaccard)'};

EmpFits = load('.\outputs\Scha400_7_lh_TopoMeasures_Thr70.mat');

sub2use = [1:298 300:973];

EmpFitType(:,1) = triu2vec(EmpFits.maxKS(sub2use,sub2use),1);
EmpFitType(:,2) = triu2vec(EmpFits.maxRMSE(sub2use,sub2use),1);
EmpFitType(:,3) = triu2vec(EmpFits.rd(sub2use,sub2use),1);
EmpFitType(:,4) = triu2vec(EmpFits.TND(sub2use,sub2use),1);
EmpFitType(:,5) = triu2vec(EmpFits.TFdiff(sub2use,sub2use),1);
EmpFitType(:,6) = triu2vec(EmpFits.DegCorr(sub2use,sub2use),1);
EmpFitType(:,7) = triu2vec(EmpFits.EdgeOverlapJacc(sub2use,sub2use),1);

save('GNN_FLaG_BestResults.mat','EmpFitType','BestFit','EmpFitType')