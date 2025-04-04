SchSize = 100:100:400;

for s = 1:length(SchSize)

Nodes = SchSize(s)/2;

verts = lh_verts;

parc_name = {['lh_scha',num2str(SchSize(s))]};

parc = [Scha7_parcs.(parc_name{1})];

Adist = zeros(Nodes);
f = waitbar(0,'Completed');
for i = 1:Nodes-1
    
    Verts_i = verts(parc==i,:);
    
    for j = i+1:Nodes
        Verts_j = verts(parc==j,:);
        
        d = triu2vec(pdist2(Verts_i,Verts_j),1);
        Adist(i,j) = mean(d);
    end
    waitbar(i/(Nodes-1),f)
end

A_dist = Adist + Adist';

save(['./data/Schaefer7_dist/Scha',num2str(SchSize(s)),'_EucDist_lh.mat'],'A_dist')

end