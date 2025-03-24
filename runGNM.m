function runGNM(MDL,ADDMULT,LAW,TIMING)
addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

A_dist = mdldata.A_dist;

adjs = mdldata.adj;

Input.NNodes = length(A_dist);

MdlPDMs = mdldata.hansen_maps;

% Convert to a Pearson distance
for i = 1:length(mdldata.hansen_maps)
    MdlPDMs{i} = 1-mdldata.hansen_maps{i};
end

nMaps = length(mdldata.hansen_maps);

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

disp(['Running model GNM_',num2str(l),'_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  

Input.Timing = l;
Output = fitGNM(adjs{1},A_dist,PDMs,Input);
Output.PDMs = PDMs;

save(['GNM_',num2str(l),'_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')
    end
end

end

end