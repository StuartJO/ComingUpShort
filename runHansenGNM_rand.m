function runHansenGNM_rand(MDL,ADDMULT,LAW,TIMING)
addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

load('Scha7_400_COGs.mat')

A_dist = mdldata.A_dist;

adjs = mdldata.adj;

Input.NNodes = length(A_dist);

load('Scha400_SArand.mat','C')
MdlPDMs{1} = 1-C;
n = 2;

Input.useParfor = 1;
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

if k == 2
    Alphas = [1 1; repmat([0 10],n-1,1)];
    Input.PDMAlphaRange = Alphas;
end

for l = TIMING

    for i = MDL
    PDMs{2} = rescale(MdlPDMs{i});

n = length(PDMs);

disp(['Running model Hansen_Timing_',num2str(l),'_RandMdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  

Input.Timing = l;
Output = GenMdl_Timing(adjs{1},A_dist,PDMs,Input);
Output.PDMs = PDMs;
save(['Hansen_Timing_',num2str(l),'_RandMdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')
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
