load('./data/surface_data/fsaverage_surface_data.mat')

verts = [lh_verts;rh_verts];

r_parc = Scha7_parcs.rh_scha400+200;

r_parc(r_parc==200) = 0;

parc = [Scha7_parcs.lh_scha400 r_parc];

A_dist_ = zeros(400);
f = waitbar(0,'Completed');
for i = 1:399
    
    Verts_i = verts(parc==i,:);
    
    for j = i+1:400
        Verts_j = verts(parc==j,:);
        
        d = pdist2(Verts_i,Verts_j);
        A_dist_(i,j) = mean(d(:));
    end
    waitbar((i+1)/400,f,[num2str((i+1)),'/400'])
end

A_dist = A_dist_+A_dist_';

save('./data/Schaefer7_dist/Scha400_EucDist_full.mat','A_dist')