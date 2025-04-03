function MdlOutput = fitGNM_FLaG(adjs,A_dist,PDMs,Input)

% This function runs the generative model optimisation. 
% The model has the following basic form for defining connection probabilities:
%
% (PD1^eta)*a1(T^gam)*a2(PD2^lam) or (PD1^eta)+a1(T^gam)+a2(PD^lam)
%
% where PD1/PD2 is a measure of distance/similarity between a pair of
% nodes; T is some measure of the topology between a pair of nodes; and
% eta, gam, lam, a1, and a2 are free parameters. Note that each power-law
% interaction (e.g., PD1^eta) can be replaced with an exponential function
% (e.g., exp(eta*PD1)).
%
% See the bottom of this header as to the model form appears different to 
% what is in the paper/why there are seemingly more than 3 free parameters
%
% Teachnically, this function is a wrapper for another function
% (GrowthModel) which itself is a wrapper for functions which actually run
% the generative model. Turtles all the way down.
%
% Inputs:
% 
% A = The adjacency matrix to model
%
% A_dist = The distances between nodes (e.g., Euclidean distances, fibre
% distance) for A. This distance is used to compute how similar the model
% is to the input adjacency matrix.
%
% PD1 = Either a matrix indicating pairwise similarity/distances between nodes, or a cell where each
% element contains such a distance matrix, to be used for the modelling itself.
% The "growth" model is specified by the latter of these options.
%
% PD2 = A second matrix indicating pairwise similarity/distances between nodes in
% A (unlike with PD1 a cell cannot be used as an input here). This could be 
% correlated gene expression, similarity in histology etc.
%
% Input = a structure containing the many possible options needed to
% configure the model. The following fields are required:
%   AddMult = 'Add' or 'Mult', if the multiplicative or additive form is to
%   be used
%
%   ModelNum = A value between 1 and 13, each corresponds to the following
%   topology form:
%       1.  'sptl'          spatial model
%       2.  'neighbors'     number of common neighbors
%       3.  'matching'      matching index
%       4.  'clu-avg'       average clustering coeff.
%       5.  'clu-min'       minimum clustering coeff.
%       6.  'clu-max'       maximum clustering coeff.
%       7.  'clu-diff'      difference in clustering coeff.
%       8.  'clu-prod'      product of clustering coeff.
%       9.  'deg-avg'       average degree
%       10. 'deg-min'       minimum degree
%       11. 'deg-max'       maximum degree
%       12. 'deg-diff'      difference in degree
%       13. 'deg-prod'      product of degree
%       If not specified, Input.TopoType will be the name of the topological form 
%
%  ParamRange = A 5x2 matrix where the first column gives the lower, and
%  the second the upper, bounds for the following parameters (see the basic
%  form outlined at the top):
%       ParamRange(1,:) = eta
%       ParamRange(2,:) = gam
%       ParamRange(3,:) = a1
%       ParamRange(4,:) = a2
%       ParamRange(5,:) = lam
%  'ParamRange(i,:) = [x x]' will mean 'x' is always used as an input for
%  parameter 'i'. 'ParamRange(i,:) = [NaN NaN]' means that parameter isn't
%  being included (if this parameter is required to be a value, you will
%  encounter an error)
%
%  pow = the severity/exponent for the opotimisation. Higher values will
%  make it more likely to sample from cells that have the smallest maxKS
%
%  nlvl = number of steps in the optimisation
%
%  ndraw = number of repetitions/samples per step
%
%  PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%  PD1 and eta. Defaults to 'exponential'
%
%  TopoFunc = 'power-law' or 'exponential'. Controls the interaction between T
%  and eta. Defaults to 'power-law'
%
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%   PD2Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD2 and lam. Defaults to 'power-law'
%
%   epsilon = an amount to add to each edges topology value in the model (to
%   ensure each edge doesn't become undefinied). Defaults to 0 when 
%   Input.AddMult = 'Add' or 1e-6 when Input.AddMult = 'Mult'
%
%   seed = a seed network. If none is desired set to [] (default)
%
%   useParfor = set to 1 to use parfor loops where possible (0 by default)
%
%   normsum = set to 1 to normalise each term by its sum (set to 0 by
%   default, which normalises by the max)
%
% Outputs:
% 
% MdlOutput = a structure with the following fields:
%   maxKS = N*1 array of maxKS values for each network generated
%       during optimisation
%   DegCorr = N*1 array of values for the correlation with  
%       empirical degree each network generated during optimisation  
%   KS = N*4 matrix of KS values for each network generated
%       during optimisation. Each column corresponds to the following 
%       measure: 1 = degree; 2 = clustering; 3 = betweenness; 4 = mean edge 
%       length
%   P = N*5 matrix of parameter values for each network generated
%       during optimisation
%   b = 1*N cell array, where each cell contains an edge index
%       list for each network generated during optimisation
% 
%   optim_maxKS = for networks generated using the parameters
%       of the lowest maxKS value, a 1*100 array of maxKS values
%   optim_KS = for networks generated using the parameters
%       of the lowest maxKS value, an 100*4 matrix of KS values
%   optim_b = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the lowest maxKS value       
%   optim_DegCorr = for networks generated using the parameters
%       of the lowest maxKS value, a 1*100 array of values for the 
%       correlation with empirical degree      
%
%   bestDegCorr_maxKS = for networks generated using the parameters
%       of the largest degree correlation, a 1*100 array of maxKS values
%   bestDegCorr_KS = for networks generated using the parameters
%       of the largest degree correlation, an 100*4 matrix of KS values
%   bestDegCorr_b = a 1*100 cell array where each cell contains an edge
%       index list for each network generated during using the parameters 
%       of the largest degree correlation    
%   bestDegCorr_DegCorr = for networks generated using the parameters
%       of the largest degree correlation, a 1*100 array of values for the 
%       correlation with empirical degree   
%   Input = the initial input settings
%
% You may wonder why the model has a different form to what we report in
% the paper e.g., (D^eta)+a(T^gam) or (D^eta)+a(PC^gam) etc. Simply,
% because a) we never run a model with distance + topology + gene, b) we
% wanted to avoid discussing too many parameters because that would get
% confusing quick and readers would find it difficult to keep track of, and
% c) topology and PC require seperate configuration of their parameters the
% way the model is coded. So while the parameters in the topological and PC
% models perform the same function theoretically, practically they are
% different. This however does mean you can run a model with all 5 
% parameters...(good luck getting the optimisation to work!)
%
% Another question you may have is why the ParamRange variable is ordered in
% this way i.e., two exponents, then the two alpha values, then an
% exponent. The simple answer is this just follows the order in which the
% function developed over time. eta and gamma are the OG parameters so they
% come first. We tried a simple model using only alpha for the topology
% term so that came third. Then we added the PD2 term, which we first tried 
% with just an alpha value (4th), then with its own exponent as well (5th).

rng('shuffle')

AddMult = Input.AddMult;
PDMfunc = Input.PDMfunc;

if ~isfield(Input,'epsilon')
    % switch Input.AddMult
    % case 'Add'
    % Input.epsilon = 0;
    % case 'Mult'
    % Input.epsilon = 1e-6;
    % end
    Input.epsilon = 1e-6;
end

if ~isfield(Input,'seed')
    Input.seed = [];  
end

if ~isfield(Input,'useParfor')
 Input.useParfor = 0;  
end

if ~isfield(Input,'GenFromBest')
 Input.GenFromBest = 0;  
end


n_subs = length(adjs);
adjs_m = zeros(n_subs,1);

for i = 1:n_subs

[~,~,adjs_m(i)] = density_und(adjs{i}); 

end

if length(unique(adjs_m)) == 1
    SpeedUp = 1;
    disp('All empirical networks have the same density :)')
else
    SpeedUp = 0;  
end

m = max(adjs_m);

n = length(A_dist);

if ~exist(Input.SaveName_temp,'file')
disp('Running model!')
[maxKS,KS,P,b,DegCorr,Input,bvals] = VoronoinLandScape_GNMFLaG(adjs,A_dist,PDMs,m,Input);

MdlOutput.maxKS = maxKS;
MdlOutput.DegCorr = DegCorr;
MdlOutput.KS = KS;
MdlOutput.P = P;
MdlOutput.b = b;
MdlOutput.bvals = bvals;

% Save the input configurations to output. Helps to keep track of what was
% done
MdlOutput.Input = Input;
MdlOutput.Input.NNodes = n;

save(Input.SaveName_temp,'-struct','MdlOutput','-v7.3')

else
disp('Load existing data!')
load(Input.SaveName_temp)

MdlOutput.maxKS = maxKS;
MdlOutput.DegCorr = DegCorr;
MdlOutput.KS = KS;
MdlOutput.P = P;
MdlOutput.b = b;

MdlOutput.Input = Input;
MdlOutput.Input.NNodes = n;
bvals = cell(length(b),1);
for i = 1:length(b)
    B = zeros(n);
    B(b{i})=1;
    B = B + B';
bvals{i}{1} = sum(B,2);
bvals{i}{2} = clustering_coef_bu(B);
bvals{i}{3} = betweenness_bin(B)';
bvals{i}{4} = A_dist(triu(B,1) > 0);

end

end

A_vals = cell(n_subs,1);
ANetStats4Eval = cell(n_subs,1);
for i = 1:n_subs
A = adjs{i};
A_vals{i}{1} = sum(A,2);
A_vals{i}{2} = clustering_coef_bu(A);
A_vals{i}{3} = betweenness_bin(A)';
A_vals{i}{4} = A_dist(triu(A,1) > 0);
ANetStats4Eval{i} = CalcNetStats4Eval(A,A_dist,A_vals{i});
end

nNets = length(b);

MdlOutput.EdgeOverlap = zeros(nNets,n_subs,4);
MdlOutput.TND = zeros(nNets,n_subs);
MdlOutput.maxRMSE = zeros(nNets,n_subs);
MdlOutput.maxRd = zeros(nNets,n_subs);
MdlOutput.TFdiff = zeros(nNets,n_subs);
MdlOutput.TFdiffnc = zeros(nNets,n_subs);

for i = 1:nNets
    B = zeros(n);
    B(b{i})=1;
    B = B + B';
    if SpeedUp
bNetStats4Eval = CalcNetStats4Eval(B,A_dist,bvals{i});
    end
for j = 1:n_subs
A = adjs{j};
    if SpeedUp == 0
        bNetStats4Eval = CalcNetStats4Eval(B,A_dist,bvals{i,j});
    end
[~,EdgeOverlap,MdlOutput.TND(i,j),MdlOutput.maxRMSE(i,j),MdlOutput.rd(i,j),MdlOutput.TFdiff(i,j),MdlOutput.TFdiffnc(i,j)] = CalcEvalStats(A,B,ANetStats4Eval{j},bNetStats4Eval,1);
MdlOutput.EdgeOverlap(i,j,:) = EdgeOverlap;
end

end

[~,I] = min(mean(maxKS,2));     
optim_P = P(I,:);

if Input.GenFromBest>0

NGenFromBest = Input.GenFromBest;

optim_b = cell(1,NGenFromBest);

optim_maxKS = zeros(NGenFromBest,n_subs);
optim_KS = zeros(NGenFromBest,4,n_subs);
optim_DegCorr = zeros(NGenFromBest,n_subs);

if Input.useParfor
    for k = 1:NGenFromBest
        [~,optim_b{k}] = GNM(PDMs,optim_P,m,Input);
        
        parfor i = 1:n_subs
            B = zeros(n);
            B(optim_b{k}(1:adjs_m(i))) = 1;
            B = B + B';
            [optim_maxKS(k,i), optim_KS(k,:,i)] = calc_maxKS(A_vals{i},A_dist,B);
            optim_DegCorr(k,i) = corr(sum(B)',A_vals{i}{1},'Type','Spearman');            
        end
    end   
else
    for k = 1:NGenFromBest
        [~,optim_b{k}] = GNM(PDMs,optim_P,m,Input);
        for i = 1:n_subs
            B = zeros(n);
            B(optim_b{k}(1:adjs_m(i))) = 1;
            B = B + B';
            [optim_maxKS(k,i), optim_KS(k,:,i)] = calc_maxKS(A_vals{i},A_dist,B);
            optim_DegCorr(k,i) = corr(sum(B)',A_vals{i}{1},'Type','Spearman');            
        end
    end     
end

[~,I] = max(mean(DegCorr,2));     
bestDegCorr_P = P(I,:);
bestDegCorr_b = cell(1,NGenFromBest);

bestDegCorr_maxKS = zeros(NGenFromBest,n_subs);
bestDegCorr_KS = zeros(NGenFromBest,4,n_subs);
bestDegCorr_DegCorr = zeros(NGenFromBest,n_subs);

if Input.useParfor
    for k = 1:NGenFromBest
        [~,bestDegCorr_b{k}] = GNM(PDMs,bestDegCorr_P,m,Input);
        
        parfor i = 1:n_subs
            B = zeros(n);
            B(bestDegCorr_b{k}(1:adjs_m(i))) = 1;
            B = B + B';
            [bestDegCorr_maxKS(k,i), bestDegCorr_KS(k,:,i)] = calc_maxKS(A_vals{i},A_dist,B);
            bestDegCorr_DegCorr(k,i) = corr(sum(B)',A_vals{i}{1},'Type','Spearman');            
        end
    end   
else
    for k = 1:NGenFromBest
        [~,bestDegCorr_b{k}] = GNM(PDMs,bestDegCorr_P,m,Input);
        for i = 1:n_subs
            B = zeros(n);
            B(bestDegCorr_b{k}(1:adjs_m(i))) = 1;
            B = B + B';
            [bestDegCorr_maxKS(k,i), bestDegCorr_KS(k,:,i)] = calc_maxKS(A_vals{i},A_dist,B);
            bestDegCorr_DegCorr(k,i) = corr(sum(B)',A_vals{i}{1},'Type','Spearman');            
        end
    end     
end
       
MdlOutput.BestFit.maxKS.maxKS = optim_maxKS;
MdlOutput.BestFit.maxKS.KS = optim_KS;
MdlOutput.BestFit.maxKS.b = optim_b;
MdlOutput.BestFit.maxKS.DegCorr = optim_DegCorr;
MdlOutput.BestFit.maxKS.P = optim_P;

MdlOutput.BestFit.Deg.maxKS = bestDegCorr_maxKS;
MdlOutput.BestFit.Deg.KS = bestDegCorr_KS;
MdlOutput.BestFit.Deg.b = bestDegCorr_b;
MdlOutput.BestFit.Deg.DegCorr = bestDegCorr_DegCorr;
MdlOutput.BestFit.Deg.P = bestDegCorr_P;

end



