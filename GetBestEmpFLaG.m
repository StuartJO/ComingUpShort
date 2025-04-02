
%Nsub = 972;
MdlBestFitIndv = zeros(10,7,7,972);
MdlBestFitAll = cell(10,7);
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
    if ismember(j,[6 7])
    [maxFit,I] = max(FIT{j},[],'linear');
    [~,BestNet] = max(FIT{j});
    else
    [minFit,I] = min(FIT{j},[],'linear');
    [~,BestNet] = min(FIT{j});
    end
    UniqueBestNet = unique(BestNet);
    MdlBestFitAll{i,j} = zeros(length(UniqueBestNet)*Nsub,7);
    for k = 1:7
        fitdata = FIT{k}(UniqueBestNet,:);
        MdlBestFitAll{i,j}(:,k) = fitdata(:);
        fitdata = FIT{k}(I);
        MdlBestFitIndv(i,j,k,:) = fitdata(:);
    end
end

end

FitName = {'max(\itKS\rm)','max(\itRMSE\rm)','max(\itr_d\rm )','\itTND','\itTF_{diff }','Degree correlation','Connection overlap (Jaccard)'};

EmpFitInput = load('.\outputs\Scha400_7_lh_TopoMeasures_Thr70.mat');

sub2use = [1:298 300:973];

EmpFit(:,1) = triu2vec(EmpFitInput.maxKS(sub2use,sub2use),1);
EmpFit(:,2) = triu2vec(EmpFitInput.maxRMSE(sub2use,sub2use),1);
EmpFit(:,3) = triu2vec(EmpFitInput.rd(sub2use,sub2use),1);
EmpFit(:,4) = triu2vec(EmpFitInput.TND(sub2use,sub2use),1);
EmpFit(:,5) = triu2vec(EmpFitInput.TFdiff(sub2use,sub2use),1);
EmpFit(:,6) = triu2vec(EmpFitInput.DegCorr(sub2use,sub2use),1);
EmpFit(:,7) = triu2vec(EmpFitInput.EdgeOverlapJacc(sub2use,sub2use),1);

Mdl_names = {'Spatial','Gene coexpression','Receptor similarity','Laminar similarity','Metabolic connectivity','Haemodynamic connectivity','Electrophysiological connectivity','Temporal similarity','Random similarity','Matching'};

save('GNM_FLaG_BestResults.mat','EmpFit','MdlBestFitIndv','FitName','Mdl_names','MdlBestFitAll','-v7.3')