
load('C:\Users\Stuart\Documents\GitHub\ComingUpShort\data\Hansen_networks.mat')

A = adj{1};

ANetStats4Eval = CalcNetStats4Eval(A,A_dist);

AddMult = {'Add','Mult'};
exponent = {'powerlaw','exponential'};

Nnodes = length(A_dist);

for i = 1:2
    addmult = AddMult{i};
    for j = 1:2     
        expo = exponent{j};
        for mdl = 0:7

Output = load(['C:\Users\Stuart\Documents\GitHub\ComingUpShort\outputs\old\Hansen_Timing_0_Mdl',num2str(mdl),'_',addmult,'_',expo,'.mat']);

FitMetricsOutput = load(['C:\Users\Stuart\Documents\GitHub\ComingUpShort\outputs\old\Hansen_Timing_0_Mdl',num2str(mdl),'_',addmult,'_',expo,'_TopoTypes.mat']); 

MdlOutput.EdgeOverlap = [FitMetricsOutput.EdgeOverlap FitMetricsOutput.EdgeJaccard];

MdlOutput.TFdiff = FitMetricsOutput.TFdiff;

MdlOutput.TFdiffnc = FitMetricsOutput.TFdiff_no_c;

MdlOutput.TND = pdist2(FitMetricsOutput.TopoGlobal,ANetStats4Eval.GlobalTopo);

MdlOutput.maxRd = 1-min(FitMetricsOutput.TopogCorr,[],2);

MdlOutput.maxRMSE = max(FitMetricsOutput.RMSE(:,1:4),[],2);

MdlOutput.maxKS = Output.maxKS;
MdlOutput.DegCorr = Output.DegCorr;
MdlOutput.KS = Output.KS;
MdlOutput.P = Output.P;
MdlOutput.b = Output.b;
MdlOutput.PDMs = Output.PDMs;
MdlOutput.Input = Output.Input;

%BestOutput = load(['C:\Users\Stuart\Documents\GitHub\ComingUpShort\outputs\old\Hansen_Timing_0_Mdl',num2str(mdl),'_',addmult,'_',expo,'_optim_FitType.mat']); 

nGenFromBest = length(Output.optim_b);

% Predefine variables with appropriate sizes:
optim_EdgeOverlap = zeros(nGenFromBest, 4);  % Adjust size if known
optim_TND = zeros(nGenFromBest, 1);         % TND likely scalar per iteration
optim_maxRMSE = zeros(nGenFromBest, 1);     % maxRMSE likely scalar per iteration
optim_rd = zeros(nGenFromBest, 1);          % rd likely scalar per iteration
optim_TFdiff = zeros(nGenFromBest, 1);      % TFdiff likely scalar per iteration
optim_TFdiffnc = zeros(nGenFromBest, 1);    % TFdiffnc likely scalar per iteration


for k = 1:nGenFromBest
    B = zeros(Nnodes);
    B(Output.optim_b{k})=1;
    B = B + B';

    bNetStats4Eval = CalcNetStats4Eval(B,A_dist);
    [~,optim_EdgeOverlap(k,:),optim_TND(k),optim_maxRMSE(k),optim_rd(k),optim_TFdiff(k),optim_TFdiffnc(k)] = CalcEvalStats(A,B,ANetStats4Eval,bNetStats4Eval,0);
end     

% Predefine variables with appropriate sizes:
bestDegCorr_EdgeOverlap = zeros(nGenFromBest, 4);  % Adjust size if known
bestDegCorr_TND = zeros(nGenFromBest, 1);         % TND likely scalar per iteration
bestDegCorr_maxRMSE = zeros(nGenFromBest, 1);     % maxRMSE likely scalar per iteration
bestDegCorr_rd = zeros(nGenFromBest, 1);          % rd likely scalar per iteration
bestDegCorr_TFdiff = zeros(nGenFromBest, 1);      % TFdiff likely scalar per iteration
bestDegCorr_TFdiffnc = zeros(nGenFromBest, 1);    % TFdiffnc likely scalar per iteration


for k = 1:100
    B = zeros(Nnodes);
    B(Output.bestDegCorr_b{k})=1;
    B = B + B';

    bNetStats4Eval = CalcNetStats4Eval(B,A_dist);

    [~,bestDegCorr_EdgeOverlap(k,:),bestDegCorr_TND(k),bestDegCorr_maxRMSE(k),bestDegCorr_rd(k),bestDegCorr_TFdiff(k),bestDegCorr_TFdiffnc(k)] = CalcEvalStats(A,B,ANetStats4Eval,bNetStats4Eval,0);

end  

      
MdlOutput.BestFit.maxKS.maxKS = Output.optim_maxKS;
MdlOutput.BestFit.maxKS.KS = Output.optim_KS;
MdlOutput.BestFit.maxKS.b = Output.optim_b;
MdlOutput.BestFit.maxKS.DegCorr = Output.optim_DegCorr;

[~,I]= min(Output.maxKS);

MdlOutput.BestFit.maxKS.P = Output.P(I,:);

MdlOutput.BestFit.maxKS.meanProb = Output.optim_meanProb;
MdlOutput.BestFit.maxKS.medianProb = Output.optim_medianProb;

MdlOutput.BestFit.DegCorr.maxKS = Output.bestDegCorr_maxKS;
MdlOutput.BestFit.DegCorr.KS = Output.bestDegCorr_KS;
MdlOutput.BestFit.DegCorr.b = Output.bestDegCorr_b;
MdlOutput.BestFit.DegCorr.DegCorr = Output.bestDegCorr_DegCorr;

[~,I]= max(Output.DegCorr);

MdlOutput.BestFit.maxKS.P = Output.P(I,:);

MdlOutput.BestFit.DegCorr.meanProb = Output.bestDegCorr_meanProb;
MdlOutput.BestFit.DegCorr.medianProb = Output.bestDegCorr_medianProb;

MdlOutput.BestFit.maxKS.EdgeOverlap = optim_EdgeOverlap;
MdlOutput.BestFit.maxKS.TND = optim_TND;
MdlOutput.BestFit.maxKS.maxRMSE = optim_maxRMSE;
MdlOutput.BestFit.maxKS.maxRd = optim_rd;
MdlOutput.BestFit.maxKS.TFdiff = optim_TFdiff;
MdlOutput.BestFit.maxKS.TFdiffnc = optim_TFdiffnc;

MdlOutput.BestFit.DegCorr.EdgeOverlap = bestDegCorr_EdgeOverlap;
MdlOutput.BestFit.DegCorr.TND = bestDegCorr_TND;
MdlOutput.BestFit.DegCorr.maxRMSE = bestDegCorr_maxRMSE;
MdlOutput.BestFit.DegCorr.maxRd = bestDegCorr_rd;
MdlOutput.BestFit.DegCorr.TFdiff = bestDegCorr_TFdiff;
MdlOutput.BestFit.DegCorr.TFdiffnc = bestDegCorr_TFdiffnc;

save(['GNM_Mdl',num2str(mdl),'_',addmult,'_',expo,'.mat'],'-struct','MdlOutput','-v7.3')

        end

    end
end

