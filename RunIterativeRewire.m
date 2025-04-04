load('Hansen_networks.mat')
A = adj{1};

RewireOrder = {'random','shortest','longest','shortest','longest'};
DistMatch = {'none','prob','prob','invprob','invprob'};

for i = 1:5

    [FitMeasures,EdgeOverlap,~] = IterativeRewire(A,A_dist,'Rewires',1973,'Repeats',30,'RewireOrder',RewireOrder{i},'DistMatchType',DistMatch{i});

end

save('Network_rewirings.mat','EdgeOverlap','FitMeasures')
