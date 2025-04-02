load('Hansen_networks.mat')
A = adj{1};

load('Hansen400_rewirings_old.mat')
GlobalA = [efficiency_bin(A) diffusion_efficiency(A) modularity_Q(A) transitivity_bu(A)];

for j = 1:5

FitMeasures{j}.maxRMSE = max(TopoMeasures{j}.RMSE(:,:,1:4),[],3);
FitMeasures{j}.TF = TopoMeasures{j}.TF;
FitMeasures{j}.DegCorr = TopoMeasures{j}.TopogCorr(:,:,1);
GlobalB = TopoMeasures{j}.GlobalTopo;
GolBSz = size(GlobalB);
% reshape(GlobalB,[prod(GolBSz(1:2)) 4]);
TopoDist_ = pdist2(GlobalA,reshape(GlobalB,[prod(GolBSz(1:2)) 4]));
FitMeasures{j}.TND = reshape(TopoDist_,GolBSz(1:2));
FitMeasures{j}.maxKS = maxKS{j};
FitMeasures{j}.maxRd = max(1-TopoMeasures{j}.TopogCorr,[],3);
FitMeasures{j}.KS = KS{j};
end

save('./data/Hansen400_rewirings.mat','FitMeasures','EdgeOverlap')
