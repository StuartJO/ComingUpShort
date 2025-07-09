

 load('.\data\compiled_outputs\BestMdls_GNM_maxKS_Add_exponential.mat')

 MmaxKS = mean(maxKS,2);
 SDmaxKS = std(maxKS,[],2);
 [~,ord] = sort(MmaxKS,'ascend');

disp('Best GNMs for maxKS')

 for i = 1:3
     I = ord(i);
disp([Mdl_names{I},': M=',num2str(MmaxKS(I)),': SD=',num2str(SDmaxKS(I))])
 end

 Movlp = mean(R(:,:,1),2);
 disp(['Overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])


  MDegCorr = mean(DegCorr,2);

 disp(['DegCorr range: ',num2str(min(MDegCorr)),' to ',num2str(max(MDegCorr))])


DistBins = {'short','mid','long'};
for i = 1:3
 Movlp = mean(R(:,:,1+i),2);
 disp([DistBins{i},' overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])
end

 load('C:\Users\Stuart\Documents\GitHub\ComingUpShort\data\compiled_outputs\BestMdls_GNM_DegCorr_Add_exponential.mat')
disp('Best GNMs for DegCorr')

 MmaxKS = mean(maxKS,2);
 SDmaxKS = std(maxKS,[],2);
 [~,ord] = sort(MmaxKS,'ascend');

 for i = 1:3
     I = ord(i);
disp([Mdl_names{I},': M=',num2str(MmaxKS(I)),': SD=',num2str(SDmaxKS(I))])
 end

 Movlp = mean(R(:,:,1),2);
 disp(['Overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])

  MDegCorr = mean(DegCorr,2);

 disp(['DegCorr range: ',num2str(min(MDegCorr)),' to ',num2str(max(MDegCorr))])

DistBins = {'short','mid','long'};
for i = 1:3
 Movlp = mean(R(:,:,1+i),2);
 disp([DistBins{i},' overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])
end

%%

disp('Best WB GNMs for maxKS')


 load('C:\Users\Stuart\Documents\GitHub\ComingUpShort\data\compiled_outputs\BestMdls_WB_GNM_maxKS_Add_exponential.mat')

 MmaxKS = mean(maxKS,2);
 SDmaxKS = std(maxKS,[],2);
 [~,ord] = sort(MmaxKS,'ascend');


 for i = 1:3
     I = ord(i);
disp([Mdl_names{I},': M=',num2str(MmaxKS(I)),': SD=',num2str(SDmaxKS(I))])
 end

  disp(['maxKS range: ',num2str(min(MmaxKS)),' to ',num2str(max(MmaxKS))])

 Movlp = mean(R(:,:,1),2);
 disp(['Overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])


  MDegCorr = mean(DegCorr,2);

 disp(['DegCorr range: ',num2str(min(MDegCorr)),' to ',num2str(max(MDegCorr))])


DistBins = {'short','mid','long'};
for i = 1:3
 Movlp = mean(R(:,:,1+i),2);
 disp([DistBins{i},' overlap range: ',num2str(min(Movlp)),' to ',num2str(max(Movlp))])
end

%%

load('.\data\Network_rewirings.mat')

Rewire_mean_maxKS = mean(FitMeasures{1}.maxKS);

load('.\data\compiled_outputs\BestMdls_GNM_maxKS_Add_exponential.mat')

matching_mean = min(mean(maxKS,2));

NRewires2MatchGNM = sum(matching_mean > Rewire_mean_maxKS);

Ovlp_NRewires2MatchGNM = EdgeOverlap{1}(1,NRewires2MatchGNM);


%% FLaG ranges

load('.\data\compiled_outputs\GNM_FLaG_BestResults.mat')

EmpFitRange = zeros(7,2);
for i = 1:7
    EmpFitRange(i,:) = [min(EmpFit(:,i)) max(EmpFit(:,i))];
end

for j = 1:7
FeatureRange = EmpFitRange(j,:);
for i = 1:10
    vec= MdlBestFitAll{i,j}(:,j);
Isin_range = vec >= FeatureRange(:,1) & vec <= FeatureRange(:,2);
InRange(i,j) = sum(Isin_range) / numel(vec) * 100;
end
end

 [~,ord] = sort(InRange(:,1),'descend');

 disp('Best maxKS overlap for FLaG results')

 for i = 1:3
     I = ord(i);
disp([Mdl_names{I},': maxKS overlap=',num2str(InRange(I,1))])
 end

disp(['FLaG: maxKS range of bottom 7=',num2str(min(InRange(ord(4:end),1))),' to ',num2str(max(InRange(ord(4:end),1)))])

disp(['FLaG: Edge Overlap range =',num2str(min(InRange(:,6))),' to ',num2str(max(InRange(:,6)))])
disp(['FLaG: DegCorr range =',num2str(min(InRange(:,7))),' to ',num2str(max(InRange(:,7)))])

%%

EmpFitInput = load('.\data\EmpiricalSimilarity\Schaefer_7net_iFOD2_acpc_lh_str70Thr_fitMetrics.mat');

sub2use = [1:298 300:973];

for i = 1:10
EmpmaxKS(i,:) = triu2vec(EmpFitInput.maxKS(sub2use,sub2use,i),1);
EmpOvlp(i,:) = triu2vec(EmpFitInput.EdgeJaccard(sub2use,sub2use,i),1);
EmpDegCorr(i,:) = triu2vec(EmpFitInput.DegCorr(sub2use,sub2use,i),1);
end

disp(['Empirical maxKS range: ',num2str(min(EmpmaxKS(:))),' to ',num2str(max(EmpmaxKS(:)))])

disp(['Empirical overlap range: ',num2str(min(EmpOvlp(:))),' to ',num2str(max(EmpOvlp(:)))])

for i = 1:10
disp(['Schaefer',num2str(i*100),' maxKS M=',num2str(mean(EmpmaxKS(i,:))),': SD=',num2str(std(EmpmaxKS(i,:)))])
end

for i = 1:10
disp(['Schaefer',num2str(i*100),' overlap M=',num2str(mean(EmpOvlp(i,:))),': SD=',num2str(std(EmpOvlp(i,:)))])
end

for i = 1:10
disp(['Schaefer',num2str(i*100),' degcorr M=',num2str(mean(EmpDegCorr(i,:))),': SD=',num2str(std(EmpDegCorr(i,:)))])
end