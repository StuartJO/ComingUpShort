
for i = 1:10
meanKS(i,1) = mean(triu2vec(squeeze(maxKS(sub2use,sub2use,i))));
meanKS(i,2) = std(triu2vec(squeeze(maxKS(sub2use,sub2use,i))));
end

for i = 1:10
Ovl(i,1) = mean(triu2vec(squeeze(data.EdgeOverlap{1}(sub2use,sub2use,i))));
Ovl(i,2) = std(triu2vec(squeeze(data.EdgeOverlap{1}(sub2use,sub2use,i))));
end

for i = 1:10
DegCorr(i,1) = mean(triu2vec(squeeze(data.DegCorr(sub2use,sub2use,i))));
DegCorr(i,2) = std(triu2vec(squeeze(data.DegCorr(sub2use,sub2use,i))));
end
