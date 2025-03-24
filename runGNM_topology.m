function runGNM_topology(MDL,ADDMULT,LAW)
addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

load('Scha7_400_COGs.mat')

A_dist = mdldata.A_dist;

adjs = mdldata.adj;

Input.NNodes = length(A_dist);

Input.useParfor = 0;
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

Input.TopoFunc = 'powerlaw';

Input.TopoParamRange = [-10 10];

disp(['Running model GNM_TopoMdl_',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  
 
Output = fitGNM(adjs{1},A_dist,PDMs,Input);
Output.PDMs = PDMs;
save(['GNM_TopoMdl_',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')


end

end

end