load('Hansen_networks.mat')
A = adj{1};

RewireOrder = {'random','shortest','longest','shortest','longest'};
DistMatch = {'none','prob','prob','invprob','invprob'};

FitMeasures = cell(1,5);
EdgeOverlap = cell(1,5);

tic
for i = 1:5
    [FitMeasures{i},EdgeOverlap{i},~] = IterativeRewire(A,A_dist,'Rewires',1973,'Repeats',30,'RewireOrder',RewireOrder{i},'DistMatchType',DistMatch{i});
    disp(i)
    toc
end

save('./data/Network_rewirings.mat','EdgeOverlap','FitMeasures')
