function [maxKS,EdgeOverlap,TND,maxRMSE,rd,TFdiff,TFdiffnc,TopoCorrs,RMSE] = CalcEvalStats(A,B,ANetStats,BNetStats,AvgMean)

%% Get maxKS

A_vals{1} = ANetStats.NodeMeasures(:,1);
A_vals{2} = ANetStats.NodeMeasures(:,2);
A_vals{3} = ANetStats.NodeMeasures(:,5);
A_vals{4} = ANetStats.EdgeDists;
B_vals{1} = BNetStats.NodeMeasures(:,1);
B_vals{2} = BNetStats.NodeMeasures(:,2);
B_vals{3} = BNetStats.NodeMeasures(:,5);
B_vals{4} = BNetStats.EdgeDists;
KS = zeros(1,4);
for j = 1:4
[~,~,KS(j)] = kstest2(A_vals{j},B_vals{j});
end
maxKS = max(KS);

%% Get EdgeOverlap

athr = triu2vec(A,1);
bthr = triu2vec(B,1);

EdgeOverlap(1,1) = sum(athr.*bthr)./sum(athr);
EdgeOverlap(1,2) = sum((athr==0).* (bthr==0))./sum(athr==0);
EdgeOverlap(1,3) = sqrt(prod(EdgeOverlap(:,1:2),2));
EdgeOverlap(1,4) = sum(athr & bthr)/sum(athr | bthr);  

%% Get TND

TND = pdist2(ANetStats.GlobalTopo,BNetStats.GlobalTopo);

%% Get RMSE
NNodes = size(ANetStats.NodeMeasures,1);

RMSE = zeros(1,5);

if AvgMean == 0

for i = 1:5
RMSE(1,i) = sqrt(sum((ANetStats.NodeMeasures(:,i)-BNetStats.NodeMeasures(:,i)).^2)/NNodes)/ANetStats.meanNodeMeasures(i);
end

else

for i = 1:5
RMSE(1,i) = sqrt(sum((ANetStats.NodeMeasures(:,i)-BNetStats.NodeMeasures(:,i)).^2)/NNodes)/mean([ANetStats.meanNodeMeasures(i) BNetStats.meanNodeMeasures(i)]);
end
    
end

maxRMSE = max(RMSE(1:4));

%% Get maxCorr

TopoCorrs = zeros(1,5);

for i = 1:5
    TopoCorrs(i) = corr(ANetStats.NodeMeasures(:,i),BNetStats.NodeMeasures(:,i));
end

rd = 1-min(TopoCorrs([1 2 4 5]));

%% Get TFdiff

TFdiff = norm(ANetStats.TF-BNetStats.TF);

ismissing = ANetStats.TF_ismissing | BNetStats.TF_ismissing;

TFdiffnc = norm(ANetStats.TF(~ismissing,~ismissing) - BNetStats.TF(~ismissing,~ismissing));




