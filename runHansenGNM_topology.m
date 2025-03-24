function runHansenGNM_topology(MDL,ADDMULT,LAW,TIMING)
addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

load('Scha7_400_COGs.mat')

A_dist = mdldata.A_dist;

adjs = mdldata.adj;

Input.NNodes = length(A_dist);

Input.useParfor = 1;
Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

PDMs{1} = A_dist;
Input.PDMfunc{1} = 'exponential';

TopoTypes = {'neighbors','matching','clu-avg','clu-min','clu-max','clu-diff','clu-prod','deg-avg','deg-min','deg-max','deg-diff','deg-prod'};

for i = MDL

    Input.TopoType = TopoTypes{i};

for k = ADDMULT
    if k == 1
    AddMult = 'Mult';
    Input.TopoAlphaRange = [1 1];
    else
    AddMult = 'Add';
    Input.TopoAlphaRange = [0 10];
    end

    Input.AddMult = AddMult;

for j = LAW
    
if j == 1
    Input.PDMfunc{1} = 'powerlaw';
    Input.PDMParamRange = [-10 0];
else
    Input.PDMfunc{1} = 'exponential';
    Input.PDMParamRange = [-2 0];
end

if k == 0
    Alphas = [1 1];
    Input.PDMAlphaRange = Alphas;
end

for l = TIMING

Input.TopoFunc = 'powerlaw';

Input.TopoParamRange = [-10 10];

n = length(PDMs);

disp(['Running model Hansen_Timing_',num2str(l),'_TopoMdl_',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  
 
Input.Timing = l;
Output = GenMdl_Timing(adjs{1},A_dist,PDMs,Input);
Output.PDMs = PDMs;
save(['Hansen_Timing_',num2str(l),'_TopoMdl_',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')


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
