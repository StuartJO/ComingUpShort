% 
subs2use = [1:298 300:973];
adj_full = zeros(432);
adj_cort = zeros(400);
adj_lh = zeros(200);
for i = 1:length(subs2use)
    A = ADJS{subs2use(i)};
    adj_full(:,:,i) = A;
    adj_cort(:,:,i) = A([1:200 217:416],[1:200 217:416]);
    adj_lh(:,:,i) = A(1:200,1:200);
end

connData = adj_lh;

Conn = mean(connData>0,3);
Thr1 = (Conn>.6);
meanSC = mean(connData,3);

meanSCThr1 = Thr1.*meanSC.*~eye(length(meanSC));

meanSCThr2 = strength_threshold(meanSCThr1,.1);


load('fsaverage_surface_data.mat')

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

deg = sum(meanSCThr2>0);

hubness = double(deg>=prctile(deg,90));

parc = Scha7_parcs.lh_scha400;

ExampleSurfacePlotFunction(surface,parc,hubness(1:200),parula(256),'Degree');

nnz(meanSCThr2)/2

deg2 = sum(adj{1});
hubness2 = double(deg2>=prctile(deg2,90)); 

ExampleSurfacePlotFunction(surface,parc,hubness2(1:200),parula(256),'Degree');
