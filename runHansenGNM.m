function runHansenGNM(MDL,ADDMULT,LAW,TIMING)
addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

load('Scha7_400_COGs.mat')

A_dist = mdldata.A_dist;

adjs = mdldata.adj;

Input.NNodes = length(A_dist);

MdlPDMs = mdldata.hansen_maps;

for i = 1:length(mdldata.hansen_maps)
    MdlPDMs{i} = 1-mdldata.hansen_maps{i};
end

nHansen = length(mdldata.hansen_maps);

for i = 1:3
coords = central_vertex_coordinates(:,i);

MdlPDMs{i+nHansen} = coords+coords';

end

for i = 1:length(mdldata.hansen_maps)
    coeff = pca(mdldata.hansen_maps{i});
    pc1 = coeff(:,1);
    MdlPDMs{i+nHansen+3}=squareform(pdist(pc1));
    pc1r = rescale(coeff(:,1));
    MdlPDMs{i+nHansen+nHansen+3} = pc1r+pc1r';
end

%deg_mid = sum(adjs{1}.*(A_dist>=30 & A_dist<90));
%f = deg_mid+deg_mid'+(10^-6);
%f = rescale(hubmid,"InputMin",0);
%F = f.*~eye(Input.NNodes);

%MdlPDMs{30} = F;

%MdlPDMs{31} = (matching(adjs{1})+(10^-6)).*~eye(Input.NNodes);


%load('.\data\Neuromaps_fsaverage164k.mat')
%load('.\fsaverage_surface_data.mat')
%NeuroMaps2Use = [28 55 72];
%for i = 1:length(NeuroMaps2Use)
%    parc_val = zeros(100,1);
%    for j = 1:200
%       parc_val(j) = nanmean(neuromap(Scha7_parcs.lh_scha400==j,NeuroMaps2Use(i))); 
%    end
%    MdlPDMs{i+31} = squareform(pdist(parc_val));
%    MdlPDMs{i+3+31} = (rescale(parc_val)+rescale(parc_val)').*~eye(Input.NNodes);
%end

n = 2;

Input.useParfor = 0;
Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

PDMs{1} = A_dist;
Input.PDMfunc{1} = 'exponential';

Input.PDMfunc{2} = 'powerlaw';

for k = ADDMULT
    if k == 1
    AddMult = 'Mult';
    else
    AddMult = 'Add';
    end

    Input.AddMult = AddMult;

for j = LAW
    
if j == 1
    Input.PDMfunc{1} = 'powerlaw';
    if k == 1
    Input.PDMParamRange = [-10 0;-10 10];
    else
    Input.PDMParamRange = [-10 0;-20 200];
    end
else
    Input.PDMfunc{1} = 'exponential';
    if k == 1
    Input.PDMParamRange = [-2 0;-10 10];
    else
    Input.PDMParamRange = [-2 0;-20 200];
    end
end

if k == 0
    Alphas = [1 1; repmat([0 10],n-1,1)];
    Input.PDMAlphaRange = Alphas;
end

for l = TIMING

    for i = MDL
    PDMs{2} = rescale(MdlPDMs{i});

n = length(PDMs);

disp(['Running model Hansen_Timing_',num2str(l),'_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  

Input.Timing = l;
Output = GenMdl_Timing(adjs{1},A_dist,PDMs,Input);
Output.PDMs = PDMs;
save(['Hansen_Timing_',num2str(l),'_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')
    end
end

end

end


%%
% clear Input
% 
% Input.useParfor = 1;
% Input.ndraw = 2000;
% Input.pow = 2;
% Input.nlvl = 5;
% Input.NNodes = length(A_dist);
% PDMs{1} = A_dist;
% 
% for k = 1:2
%     if k == 1
%     AddMult = 'Mult';
%     else
%     AddMult = 'Add';
%     end
% 
%     Input.AddMult = AddMult;
% 
% for j = 1:2
%     
% if j == 1
%     Input.PDMfunc{1} = 'powerlaw';
%     Input.ParamRange = [-10 0];
% else
%     Input.PDMfunc{1} = 'exponential';
%     Input.ParamRange = [-2 0];
% end
% 
% for l = 0:1
% 
% n = length(PDMs);
% 
% disp('Running Sptl model')  
% 
% Input.Timing = l;
% Output = GenMdl_Timing(adjs{1},A_dist,PDMs,Input);
% Output.PDMs = PDMs;
% save(['Hansen_Timing_',num2str(l),'_Mdl',num2str(0),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')
% 
% end
% 
% end
% 
% end
