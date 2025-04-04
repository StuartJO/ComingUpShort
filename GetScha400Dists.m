load('fsaverage_surface_data.mat')

verts = [lh_verts;rh_verts];

r_parc = Scha7_parcs.rh_scha400+200;

r_parc(r_parc==100) = 0;

parc = [Scha7_parcs.lh_scha400];

Adist_inc0 = zeros(201);
f = waitbar(0,'Completed');
for i = 0:199
    
    Verts_i = verts(parc==i,:);
    
    for j = i+1:200
        Verts_j = verts(parc==j,:);
        
        d = triu2vec(pdist2(Verts_i,Verts_j),1);
        Adist_inc0(i+1,j+1) = mean(d);
    end
    waitbar((i+1)/200,f,[num2str((i+1)),'/200'])
end

Adist_inc0_ = Adist_inc0+Adist_inc0';

HemiDist = Adist_inc0(1,2:end)+Adist_inc0(1,2:end)';

A_dist = [Adist_inc0_(2:end,2:end) HemiDist;HemiDist Adist_inc0_(2:end,2:end)];

save('Scha400_EucDist_full.mat','A_dist')