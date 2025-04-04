SchSize = 100:100:1000;
n_sub=973;
run = 1;
if run

    for x = 1:2
    
maxKS = zeros(n_sub,n_sub,length(SchSize));
DegCorr = zeros(n_sub,n_sub,length(SchSize));
TND = zeros(n_sub,n_sub,length(SchSize));
maxRMSE = zeros(n_sub,n_sub,length(SchSize));
maxRd = zeros(n_sub,n_sub,length(SchSize));
EdgeJaccard = zeros(n_sub,n_sub,length(SchSize));
TFdiffnc = zeros(n_sub,n_sub,length(SchSize));
TFdiff = zeros(n_sub,n_sub,length(SchSize));

NetStats4Eval = cell(n_sub,1);
A = cell(n_sub,1);

for s = 1:length(SchSize)
    tic
if x == 2
load(['Schaefer',num2str(SchSize(s)),'_7net_FACT_acpc_lh_wei.mat'])
else
load(['Schaefer',num2str(SchSize(s)),'_7net_iFOD2_acpc_lh_wei.mat'])    
end
load(['Scha',num2str(SchSize(s)),'_EucDist_lh.mat'])

N = SchSize(s)/2;

for i = 1:n_sub; den(i) = density_und(ADJS{i});end
den(299)=NaN;
denThr = nanmin(den)*.7;

for i = 1:n_sub

A{i} = double(strength_threshold(ADJS{i},denThr)>0);

NetStats4Eval{i} = CalcNetStats4Eval(A{i},A_dist);

end

for i = 1:n_sub-1
    for j = i+1:n_sub
        
 [maxKS(i,j,s),EdgeOverlap,TND(i,j,s),maxRMSE(i,j,s),maxRd(i,j,s),TFdiff(i,j,s),TFdiffnc(i,j,s),TopoCorrs] = CalcEvalStats(A{i},A{j},NetStats4Eval{i},NetStats4Eval{j},1);              
       
    maxKS(j,i,s)=maxKS(i,j,s);
    DegCorr(j,i,s)= corr(NetStats4Eval{i}.NodeMeasures(:,1),NetStats4Eval{j}.NodeMeasures(:,1),'Type','Spearman');
    DegCorr(i,j,s)=DegCorr(j,i,s);
    TFdiff(j,i,s) = TFdiff(i,j,s);
    TFdiffnc(j,i,s) = TFdiffnc(i,j,s);
    TND(j,i,s) = TND(i,j,s);
    maxRd(j,i,s) = maxRd(i,j,s);
    EdgeJaccard(j,i,s) = EdgeOverlap(4);
    EdgeJaccard(i,j,s) = EdgeOverlap(4);
    maxRMSE(j,i,s) = maxRMSE(i,j,s);
    end
end

disp(num2str(s))
toc
end
if x == 2
    save('Schaefer_7net_FACT_acpc_lh_str70Thr_fitMetrics.mat','maxKS','DegCorr','TFdiffnc','TFdiff','maxRd','maxRMSE','TND','EdgeJaccard')
else
    save('Schaefer_7net_iFOD2_acpc_lh_str70Thr_fitMetrics.mat','maxKS','DegCorr','TFdiffnc','TFdiff','maxRd','maxRMSE','TND','EdgeJaccard')
end
    end
end