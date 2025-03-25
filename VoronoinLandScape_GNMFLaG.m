function [maxKS,KS,P,b,DegCorr,Input,bvals] = VoronoinLandScape_GNMFLaG(adjs,A_dist,PDMs,m,Input)

% This function will perform optimisation using a Voronoi tessellation
% approach, where from an initial set of random points in parameter space,
% Voronoi tessellation will be performed to divide up the parameter space
% into cells. Each cell has an associated maxKS value (based on the
% parameter value used to draw that cell). The algorithm proceeds by
% preferentially sampling parameter values from cells with the smallest 
% maxKS, and then performing Voronoi tessellation again. This repeats a
% number of times.
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
% The "growth" model is specied by the latter of these options.
%
% PD2 = A second matrix indicating pairwise similarity/distances between nodes in
% A (unlike with PD1 a cell cannot be used as an input here). This could be 
% correlated gene expression, similarity in histology etc.
%
% m = the number of edges for the model to form. When doing a growth model
% this should be a vector where each elemenet specifies the number of edges
% which should be present in the network at the end of that timestep (it is
% NOT the number of edges to form at each step, but rather the cumulative
% edges). Can also be specified as the density value
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
%
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
%   The following are optional specifications for Input (will be set to
%   defaults if not set):
%
%   PD1Func = 'power-law' or 'exponential'. Controls the interaction between
%   PD1 and eta. Defaults to 'exponential'
%
%   TopoFunc = 'power-law' or 'exponential'. Controls the interaction between T
%   and eta. Defaults to 'power-law'
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

AddMult = Input.AddMult;

PDMfunc = Input.PDMfunc;

Timing = Input.Timing;

if ~isfield(Input,'epsilon')
    % switch Input.AddMult
    % case 'Add'
    % Input.epsilon = 0;
    % case 'Mult'
    % Input.epsilon = 1e-6;
    % end
    Input.epsilon = 1e-6;
end

if ~isfield(Input,'TopoType')
	Input.TopoType = 'none';
else
    if ~isfield(Input,'TopoFunc')
       error('A function type needs to be specified for topology') 
    end
end

if strcmp(Input.TopoType,'none')
	Input.TopoFunc = 'none';
end

if ~isfield(Input,'seed')
    Input.seed = [];  
end

if ~isfield(Input,'useParfor')
 Input.useParfor = 0;  
end

numProperties = 4;

pow = Input.pow;
nlvl = Input.nlvl;
ndraw = Input.ndraw;

useParfor = Input.useParfor;

totalsamples = ndraw*nlvl;

if ~isfield(Input,'PDMParamRange')
	Input.PDMParamRange = [];
    ParamType = [];
else
   ParamType = ones(1,size(Input.PDMParamRange,1));
end

if ~isfield(Input,'PDMAlphaRange')
	Input.PDMAlphaRange = [];
else
   ParamType = [ParamType ones(1,size(Input.PDMAlphaRange,1))*2];
end

if ~isfield(Input,'TopoParamRange')
	Input.TopoParamRange = [];
else
   ParamType = [ParamType ones(1,size(Input.TopoParamRange,1))*3];
end

if ~isfield(Input,'TopoAlphaRange')
	Input.TopoAlphaRange = [];
else
   ParamType = [ParamType ones(1,size(Input.TopoAlphaRange,1))*4];
end

Input.ParamType = ParamType;

Input.ParamRange = [Input.PDMParamRange;Input.PDMAlphaRange;...
    Input.TopoParamRange;Input.TopoAlphaRange];

bounds = Input.ParamRange;
nParams = size(bounds,1);

n_sub = length(adjs);

% Initalise output variables
P = zeros(totalsamples,nParams);
maxKS = zeros(totalsamples,n_sub);
DegCorr = zeros(totalsamples,n_sub);
KS = zeros(totalsamples,numProperties,n_sub);
powvals = linspace(0,pow,nlvl);
b = cell(1,totalsamples);


% Pull out the different parameters and sample at random from their ranges
% to get a set of points to start the optimisation from

for i = 1:nParams
    P(1:ndraw,i) = unifrnd(bounds(i,1),bounds(i,2),ndraw,1);
end

A_vals = cell(n_sub,1);

n = length(A_dist);

ADeg = zeros(n_sub,n);

adjs_m = zeros(n_sub,1);

for i = 1:n_sub
A = double(adjs{i}>0);
A_vals{i}{1} = sum(A,2);
A_vals{i}{2} = clustering_coef_bu(A);
A_vals{i}{3} = betweenness_bin(A)';
A_vals{i}{4} = A_dist(triu(A,1) > 0);

ADeg(i,:) = sum(A);
[~,~,adjs_m(i)] = density_und(adjs{i}); 
end

if length(unique(adjs_m)) == 1
    SpeedUp = 1;
    disp('All empirical networks have the same density :)')
else
    SpeedUp = 0;  
end
bvals = cell(totalsamples,n_sub);

for ilvl = 1:nlvl
    fprintf('level %i of %i\n',ilvl,nlvl);
    tic
    pow = powvals(ilvl);

    if ilvl==1
        ptsnew = P(1:ndraw,:);
    else
        ind = 1:ndraw*(ilvl - 1);
        mean_maxKS = mean(maxKS,2);

        if nParams == 1
            fcn_voronoi_select_sptl(P(ind,1),mean_maxKS(ind),ndraw,bounds,pow);
        elseif nParams <= 4
            ptsnew = fcn_voronoi_selectn(P(ind,:),mean_maxKS(ind),ndraw,bounds,pow);
        else
            ptsnew = fcn_custom_select(P(ind,:),mean_maxKS(ind),ndraw,pow);
        end
    end

    indnew = (1 + (ilvl - 1)*ndraw):(ilvl*ndraw);
    
    if useParfor
       
        maxKSpar = zeros(ndraw,n_sub);
        KSpar = zeros(ndraw,4,n_sub);
        Cpar = zeros(ndraw,n_sub);
        btemp = cell(1,ndraw);
        B_vals_temp = cell(ndraw,n_sub);

        parfor i = 1:ndraw          
                [~, btemp{i}] = GNM(PDMs,ptsnew(i,:),m,Input); 
                  
                if SpeedUp
                b_vals = cell(1,4); 
                B = zeros(n);
                B(btemp{i}) = 1;
                B = B + B';
                b_vals{1} = sum(B,2);
                b_vals{2} = clustering_coef_bu(B);
                b_vals{3} = betweenness_bin(B)';
                b_vals{4} = A_dist(triu(B,1) > 0);

                for j = 1:n_sub
                
                    [maxKSpar(i,j),KSpar(i,:,j)] = calc_maxKS(A_vals{j},A_dist,b_vals);
                    Cpar(i,j) = corr(sum(B)',ADeg(j,:)','Type','Spearman');
                    % So the parfor loop stops complaining
                    B_vals_temp{i,j} = b_vals;

                end
                
                else

                for j = 1:n_sub
                B = zeros(n);
                B(btemp{i}(1:adjs_m(j))) = 1;
                B = B + B';
                [maxKSpar(i,j),KSpar(i,:,j),~,~,~,B_vals_temp{i,j}] = calc_maxKS(A_vals{j},A_dist,B);
                Cpar(i,j) = corr(sum(B)',ADeg(j,:)','Type','Spearman');
                end

                end
        end
        
        maxKS(indnew,:) = maxKSpar;
        KS(indnew,:,:) = KSpar;
        DegCorr(indnew,:) = Cpar;
        b(indnew) = btemp;
        bvals(indnew,:) = B_vals_temp;
    else

        % for i = 1:ndraw
        %     indHere = indnew(i);
        %     [B, b{indHere}] = GenMdlVaryParam(PDMs,ptsnew(i,:),m,Input) ;
        %     [maxKS(indHere),KS(indHere,:)] = calc_maxKS(A_vals,A_dist,B);
        %     DegCorr(indHere) = corr(sum(B)',ADeg','Type','Spearman');
        % end
        maxKSpar = zeros(ndraw,n_sub);
        KSpar = zeros(ndraw,4,n_sub);
        Cpar = zeros(ndraw,n_sub);
        btemp = cell(1,ndraw);
        B_vals_temp = cell(ndraw,n_sub);

        for i = 1:ndraw          
                [~, btemp{i}] = GNM(PDMs,ptsnew(i,:),m,Input); 
                  
                if SpeedUp
                b_vals = cell(1,4); 
                B = zeros(n);
                B(btemp{i}) = 1;
                B = B + B';
                b_vals{1} = sum(B,2);
                b_vals{2} = clustering_coef_bu(B);
                b_vals{3} = betweenness_bin(B)';
                b_vals{4} = A_dist(triu(B,1) > 0);
                    for j = 1:n_sub
                    
                    [maxKSpar(i,j),KSpar(i,:,j)] = calc_maxKS(A_vals{j},A_dist,b_vals);
                    Cpar(i,j) = corr(sum(B)',ADeg(j,:)','Type','Spearman');
                    B_vals_temp{i,j} = b_vals;
                    end
                    
                else

                    for j = 1:n_sub
                    B = zeros(n);
                    B(btemp{i}(1:adjs_m(j))) = 1;
                    B = B + B';
                    [maxKSpar(i,j),KSpar(i,:,j),~,~,~,B_vals_temp{i,j}] = calc_maxKS(A_vals{j},A_dist,B);
                    Cpar(i,j) = corr(sum(B)',ADeg(j,:)','Type','Spearman');
                    end

                end
        end
        
        maxKS(indnew,:) = maxKSpar;
        KS(indnew,:,:) = KSpar;
        DegCorr(indnew,:) = Cpar;
        b(indnew) = btemp;
        bvals(indnew,:) = B_vals_temp;
    end
    P(indnew,:) = ptsnew;
    toc
end

if SpeedUp
    bvals = bvals(:,1);
end