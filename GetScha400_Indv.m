load(['Schaefer400_7net_iFOD2_acpc_lh_wei.mat'])    

ADJS(299)=[];
SUBS(299)=[];
n_sub = length(ADJS);
for i = 1:n_sub
    den(i) = density_und(ADJS{i});
end

denThr = nanmin(den)*.7;

for i = 1:n_sub
    THR{i} = double(strength_threshold(ADJS{i},denThr)>0);
end

save('./data/Indv_Scha400_7_lh_Thr70.mat','THR','SUBS')
