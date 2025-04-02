load('Hansen_networks.mat')
A = adj{1};

RewireOrder = {'random','shortest','longest','shortest','longest'};
DistMatch = {'none','prob','prob','invprob','invprob'};

for i = 1:5

    [FitMeasures,EdgeOverlap,DistRewires] = IterativeRewire(A,A_dist,'Rewires',1973,'Repeats',30,'RewireOrder',RewireOrder{i},'DistMatchType',DistMatch{i});


[maxKS{i},EdgeOverlap{i},KS{i},DistRewires{i},RewiredDists{i},TopoMeasures{i}] = KSrewiring(A,A_dist,'Rewires',1973,'Repeats',30,'RewireOrder',RewireOrder{i},'DistMatchType',DistMatch{i});

end

save('Hansen400_rewirings.mat','EdgeOverlap','maxKS','KS','TopoMeasures')
