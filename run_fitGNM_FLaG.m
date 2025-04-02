function run_fitGNM_FLaG(MDL,ADDMULT,LAW)

addpath(genpath('./'))

mdldata = load('Hansen_networks.mat');

load('Scha400_7_lh_Thr70.mat')

A_dist = mdldata.A_dist;

MdlPDMs = mdldata.hansen_maps;

for i = 1:length(mdldata.hansen_maps)
    MdlPDMs{i} = 1-mdldata.hansen_maps{i};
end

nHansen = length(mdldata.hansen_maps);
load('Scha400_SArand.mat','C');
MdlPDMs{nHansen+1} = 1-C;

Input.NNodes = length(A_dist);
Input.useParfor = 0;
Input.ndraw = 2000;
Input.pow = 2;
Input.nlvl = 5;

Input.normType = 'max';

for k = ADDMULT
    if k == 1
    AddMult = 'Mult';
    else
    AddMult = 'Add';
    end

    Input.AddMult = AddMult;

    for i = MDL
    clear PDMs
    
    PDMs{1} = A_dist;
    if ismember(i,2:9)    

    PDMs{2} = rescale(MdlPDMs{i-1});
    end
    n = length(PDMs);

for j = LAW

if i == 1

    if j == 1
        Input.PDMfunc{1} = 'powerlaw';
        Input.PDMParamRange = [-10 0];
    else
        Input.PDMfunc{1} = 'exponential';
        Input.PDMParamRange = [-2 0];
    end

    Input.TopoType = 'none';
    Input.TopoFunc = 'powerlaw';

elseif ismember(i,2:9)
    Input.PDMfunc{2} = 'powerlaw';
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

    Input.TopoType = 'none';
    Input.TopoFunc = 'powerlaw';
    
    if k ~= 1
        Alphas = [1 1; repmat([0 10],n-1,1)];
        Input.PDMAlphaRange = Alphas;
    end

elseif i == 10

    if j == 1
        Input.PDMfunc{1} = 'powerlaw';
        Input.PDMParamRange = [-10 0];
    else
        Input.PDMfunc{1} = 'exponential';
        Input.PDMParamRange = [-2 0];
    end

    Input.TopoType = 'matching';
    Input.TopoFunc = 'powerlaw';
    Input.TopoParamRange = [-10 10];
   if k == 1
    AddMult = 'Mult';
    Input.TopoAlphaRange = [1 1];
    else
    AddMult = 'Add';
    Input.TopoAlphaRange = [0 10];
    end

end
    
    disp(['Running model EmpFLaG_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1}])  
    Input.SaveName_temp = ['EmpFLaG_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'_temp.mat'];

    Output = GenMdlFast_FLaG2(THR,A_dist,PDMs,Input);
    Output.PDMs = PDMs;
    save(['EmpFLaG_Mdl',num2str(i),'_',AddMult,'_',Input.PDMfunc{1},'.mat'],'-struct','Output','-v7.3')

    Input.PDMfunc = [];

end

end

end