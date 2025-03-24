

SC = readmatrix('consensusSC_wei.csv','NumHeaderLines',1);

MAP_names = {'gene_coexpression','receptor_similarity','laminar_similarity','metabolic_connectivity','haemodynamic_connectivity','electrophysiological_connectivity','temporal_similarity'};

for i = 1:7
MAPS_full{i} = readmatrix([MAP_names{i},'.csv'],'NumHeaderLines',1);
hansen_maps{i} =  MAPS_full{i}(1:200,1:200);
end

adj{1} = double(SC(1:200,1:200)>0);

load('Scha400_EucDist_full.mat')

A_dist = A_dist(1:200,1:200);

save('Hansen_networks.mat','adj','hansen_maps','A_dist');





% surf_dist = readmatrix('surface_distance.csv','NumHeaderLines',1);
% 
% SC_L = SC(1:200,1:200);
% CGE_L = CGE(1:200,1:200);
% 
% scatter(SC_L(:),CGE_L(:))
% 
% deg = sum(SC>0);
% 
% surface.vertices = lh_inflated_verts;
% surface.faces = lh_faces;
% ExampleSurfacePlotFunction(surface,Scha7_parcs.lh_scha400,deg(1:200),parula(256),'Degree');
% 
