SC = readmatrix('./data/Schaefer_hansen/consensusSC_wei.csv','NumHeaderLines',1);

MAP_names = {'gene_coexpression','receptor_similarity','laminar_similarity','metabolic_connectivity','haemodynamic_connectivity','electrophysiological_connectivity','temporal_similarity'};

for i = 1:7
MAPS_full{i} = readmatrix([MAP_names{i},'.csv'],'NumHeaderLines',1);
hansen_maps{i} =  MAPS_full{i}(1:200,1:200);
end

adj{1} = double(SC(1:200,1:200)>0);

load('./data/Schaefer7_dist/Scha400_EucDist_full.mat','A_dist')

A_dist = A_dist(1:200,1:200);

save('Hansen_networks.mat','adj','hansen_maps','A_dist');

for i = 1:7
MAPS_full{i} = readmatrix([MAP_names{i},'.csv'],'NumHeaderLines',1);
hansen_maps{i} =  MAPS_full{i};
end

adj{1} = double(SC>0);

load('./data/Schaefer7_dist/Scha400_EucDist_full.mat','A_dist')

save('Hansen_networks_WB.mat','adj','hansen_maps','A_dist');
