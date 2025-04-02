function [adj, b, Edgeprob] = GNM(PDMs,params,m,Input) 

% GNM - run a generative network model based on probabilistic wiring rules
%
% Syntax:
%   [adj, b, Edgeprob] = GNM(PDMs, params, m, Input)
%
% Inputs:
%   PDMs   - Either a cell array of pairwise distance matrices (PDMs) or a 3D array, 
%            where each PDM determines connection probability contributions.
%   params - Parameter vector containing scaling and exponent parameters for PDMs 
%            and topological terms, structured according to Input.ParamType.
%   m      - Total number of edges to add to the network.
%   Input  - Structure containing:
%            - PDMfunc: Cell array of strings ('powerlaw' or 'exponential') for each PDM.
%            - normType: Normalization type ('none', 'max', or 'sum').
%            - AddMult: Combination rule for PDM and topology probabilities ('Add' or 'Mult').
%            - epsilon: Small constant to avoid numerical issues.
%            - TopoType: Topological constraint model type (e.g., 'none', 'clu-avg', etc.).
%            - TopoFunc: Topological function form ('powerlaw' or 'exponential').
%            - ParamType: Vector indicating which elements in params correspond to PDM 
%                        scaling, PDM exponents, topology scaling, and topology exponents.
%                           
%            - NNodes: Number of nodes in the network.
%
% Outputs:
%   adj      - Symmetric adjacency matrix of the generated network.
%   b        - Vector of linear indices corresponding to the added edges.
%   Edgeprob - Matrix of edge probabilities at each step of edge formation.
%
% Description:
%   This function generates a probabilistic network by sequentially adding edges
%   according to the combined influence of multiple distance-based probability matrices 
%   (PDMs) and optional topological reinforcement. Probabilities can be combined additively 
%   or multiplicatively, and are normalized according to user specifications. 

PDMfunc = Input.PDMfunc;
% VaryEtaWTime = Input.VaryEtaWTime;
normType = Input.normType;
AddMult = Input.AddMult;
epsilon = Input.epsilon;
TopoType = Input.TopoType;
TopoFunc = Input.TopoFunc;

if iscell(PDMs)
    nPDMs = length(PDMs);
    n = length(PDMs{1});   
else
    nPDMs = size(PDMs,3);
    n = size(PDMs,1);
    pdms = PDMs;
    PDMs = cell(nPDMs,1);
    for i = 1:nPDMs
       PDMs{i} = pdms(:,:,i);
    end
end

switch TopoType
    case 'none'
        DoTopo = 0;
    otherwise
        DoTopo = 1;
end

PDMParamInd = find(Input.ParamType==1);

PDMparams = params(PDMParamInd);

if isempty(PDMparams) 
    PDMparams = ones(nPDMs,1);
end

PDMAlphaInd = find(Input.ParamType==2);

PDMAlphas = params(PDMAlphaInd);

if isempty(PDMAlphas) 
    PDMAlphas = ones(nPDMs,1);
end

TopoParamInd = find(Input.ParamType==3);

TopoParams = params(TopoParamInd);

if isempty(TopoParams) 
    TopoParams = 1;
end

TopoAlphaInd = find(Input.ParamType==4);

TopoAlpha = params(TopoAlphaInd);

if isempty(TopoAlpha) 
    TopoAlpha = 1;
end

mseed = 0;
[u,v] = find(triu(ones(n),1));
indx = ((v - 1)*n + u)';

% eta = params(1,1);
%
% if VaryEtaWTime == 1
% 
% TimeSteps = length((mseed + 1):m);
% 
% eta_vals = linspace(0,eta,TimeSteps);
% 
% PDMparams(1) = eta_vals(1);
% 
% elseif VaryEtaWTime == 2
% 
% TimeSteps = length((mseed + 1):m);
% 
% eta_vals = linspace(-eta,eta,TimeSteps);
% 
% PDMparams(1) = eta_vals(1);
% 
% end

% The code only operates on the upper triangle, so all the edges are stored
% as a vector
a = ones(1,length(indx));

% Predefine the index of edges to add
b = zeros(m,1);

nPossibleEdges = length(indx);

nPDMParams = length(PDMs);

d = zeros(nPDMParams,nPossibleEdges);
Df = zeros(nPDMParams,nPossibleEdges);
PDf = zeros(nPDMParams,nPossibleEdges);

for i = 1:nPDMParams
   d(i,:) = PDMs{i}(indx); 
end

for j = 1:nPDMParams
    mv1 = PDMfunc{j};
       
    switch mv1
        case 'powerlaw'
            Df(j,:) = d(j,:).^PDMparams(j);
        case 'exponential'       
            Df(j,:) = exp(PDMparams(j).*(d(j,:)));    
    end
    switch normType
        case 'none'
            PDf(j,:) = PDMAlphas(j).*Df(j,:);
        case 'max'
            %PDf(j,:) = PDMAlphas(j).*(Df(j,:)./max(Df(j,:).*a,[],'all'));
            PDf(j,:) = PDMAlphas(j).*(Df(j,:)./max(Df(j,:).*a));
        case 'sum'
            PDf(j,:) = PDMAlphas(j).*(Df(j,:)./sum(Df(j,:).*a,'all'));        
    end
end

if DoTopo
        % The topological functions require an adjacency matrix
        A = zeros(Input.NNodes);
        TopoData = GetKSeed(A,TopoType);         
        TopoData.u = u;
        TopoData.v = v;
         switch TopoFunc
        case 'powerlaw'
            TopoData.Fk = (TopoData.K+epsilon).^TopoParams;
        case 'exponential'       
            TopoData.Fk = exp(TopoParams.*(TopoData.K+epsilon));    
        end
         Fk = TopoData.Fk(indx);        
         switch normType
            case 'none'  
                FK = TopoAlpha*Fk;
            case 'max'
                %FK = TopoAlpha*(Fk./max(Fk.*a,[],'all'));
                FK = TopoAlpha*(Fk./max(Fk.*a));
            case 'sum'
                FK = TopoAlpha*(Fk./sum(Fk.*a,'all'));        
         end

         switch AddMult
            case 'Add'
                prob = sum([PDf; FK],1);
            case 'Mult'
                prob = prod([PDf; FK],1);
        end
         
     else
         switch AddMult
        case 'Add'
            prob = sum(PDf,1);
        case 'Mult'
            prob = prod(PDf,1);
        end
end

P = prob.*a;

EdgeInds = (mseed + 1):m;

Nsteps = length((mseed + 1):m);

Edgeprob = zeros(Nsteps,nPossibleEdges);

for i = 1:Nsteps

EdgeInd = EdgeInds(i);

Edgeprob(i,:) = P;

    C = [0 cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(EdgeInd) = r;
    
    if r<= 0
       error('The code wants to form an edge that is not possible, probably because a NaN or inf has been calculated somewhere') 
       %save('DEBUG.mat') 
    end
    a(r) = 0;
      % if VaryEtaWTime > 0  
      %   switch PDMfunc{1}
      %       case 'powerlaw'
      %           Df(1,:) = d(1,:).^eta_vals(i);
      %       case 'exponential'       
      %           Df(1,:) = exp(eta_vals(i)*(d(1,:)));    
      %   end
      % end
     for j = 1:nPDMParams
        switch normType
        case 'none'  
            PDf(j,:) = PDMAlphas(j)*Df(j,:);
        case 'max'
            %PDf(j,:) = PDMAlphas(j)*(Df(j,:)./max(Df(j,:).*a,[],'all'));
            PDf(j,:) = PDMAlphas(j)*(Df(j,:)./max(Df(j,:).*a));
        case 'sum'
            PDf(j,:) = PDMAlphas(j)*(Df(j,:)./sum(Df(j,:).*a,'all'));        
        end
     end
    
     if DoTopo
         
         TopoData = CalcTopoForm(TopoType,TopoFunc,r,TopoData,TopoParams,epsilon);
         Fk = TopoData.Fk(indx);
         
         switch normType
            case 'none'  
                FK = TopoAlpha*Fk;
            case 'max'
                %FK = TopoAlpha*(Fk./max(Fk.*a,[],'all'));
                maxFK = max(Fk.*a);
                if maxFK==0
                    FK = TopoAlpha*Fk;
                else
                FK = TopoAlpha*(Fk./maxFK);
                end
            case 'sum'
                sumFK = sum(Fk.*a,'all');
                if sumFK == 0
                FK = TopoAlpha*Fk;
                else
                FK = TopoAlpha*(Fk./sumFK);     
                end
         end
         
         switch AddMult
            case 'Add'
                prob = sum([PDf; FK],1);
            case 'Mult'
                prob = prod([PDf; FK],1);
        end
         
     else
         switch AddMult
        case 'Add'
            prob = sum(PDf,1);
        case 'Mult'
            prob = prod(PDf,1);
        end
     end
     
P = prob.*a;

end

b = indx(b);

adj = zeros(n);
adj(b) = 1;
adj = adj + adj';

function TopoData = CalcTopoForm(modeltype,modelform,r,TopoData,gam,epsilon)

% CalcTopoForm - Compute topological coupling probabilities

% Note that all this code is designed so that it only updates the values
% for edges that need updating (saves on computation time, I think...).
% This is why some of the calculations look a bit more complicated than you expect 

switch modeltype
	case {'clu-avg','clu-diff','clu-max','clu-min','clu-prod'}

	c = TopoData.c;
	K = TopoData.K;
    k = TopoData.k;
	Fk = TopoData.Fk;
	A = TopoData.A;
	u = TopoData.u;
	v = TopoData.v;

	uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    k([uu,vv]) = k([uu,vv]) + 1;
    bu = A(uu,:);
    su = A(bu,bu);
    bv = A(vv,:);
    sv = A(bv,bv);
    bth = bu & bv;
    c(bth) = c(bth) + 2./(k(bth).^2 - k(bth));
    c(uu) = nnz(su)/(k(uu)*(k(uu) - 1));
    c(vv) = nnz(sv)/(k(vv)*(k(vv) - 1));
    c(k <= 1) = 0;
    bth([uu,vv]) = true;
    switch modeltype
        case 'clu-avg'
    
            K(:,bth) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')/2 + epsilon;
            K(bth,:) = bsxfun(@plus,c(:,ones(1,sum(bth))),c(bth,:)')'/2 + epsilon;

        case 'clu-diff'

            K(:,bth) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)')) + epsilon;
            K(bth,:) = abs(bsxfun(@minus,c(:,ones(1,sum(bth))),c(bth,:)'))' + epsilon;

        case 'clu-max'
            
            K(:,bth) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@max,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
            
        case 'clu-min'
            
            K(:,bth) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)') + epsilon;
            K(bth,:) = bsxfun(@min,c(:,ones(1,sum(bth))),c(bth,:)')' + epsilon;
    
        case 'clu-prod'
            
            K(bth,:) = (c(bth,:)*c') + epsilon;
            K(:,bth) = (c*c(bth,:)') + epsilon;
    end

    switch modelform
        case 'powerlaw'
            Fk(bth,:) = ((K(bth,:)).^gam);
            Fk(:,bth) = ((K(:,bth)).^gam);
        case 'exponential'
            Fk(bth,:) = exp((K(bth,:))*gam);
            Fk(:,bth) = exp((K(:,bth))*gam);
    end

    TopoData.k = k;
	TopoData.Fk = Fk;
	TopoData.c = c;

    %%

    case {'deg-avg','deg-diff','deg-max','deg-min','deg-prod'}

	k = TopoData.k;
	Fk = TopoData.Fk;
	u = TopoData.u;
	v = TopoData.v;

    w = [u(r),v(r)];
    k(w) = k(w) + 1;
    
    switch modeltype
        
    case 'deg-avg'
        
        switch modelform
        case 'powerlaw'
            Fk(:,w) = [((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam;
            Fk(w,:) = ([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam);
            Fk(w,:) = exp([((k + k(w(1)))/2) + epsilon, ((k + k(w(2)))/2) + epsilon]*gam)';
        end    

    case 'deg-diff'
        switch modelform
        case 'powerlaw'
            Fk(:,w) = (abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam;
            Fk(w,:) = ((abs([k - k(w(1)), k - k(w(2))]) + epsilon).^gam)';
        case 'exponential'
            Fk(:,w) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam);
            Fk(w,:) = exp((abs([k - k(w(1)), k - k(w(2))]) + epsilon)*gam)';
        end

    case 'deg-min'

        switch modelform
        case 'powerlaw'
            Fk(:,w) = [min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([min(k,k(w(1))) + epsilon, min(k,k(w(2))) + epsilon]*gam)';
        end

    case 'deg-max'
        
        switch modelform
        case 'powerlaw'
            Fk(:,w) = [max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam;
            Fk(w,:) = ([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon].^gam)';
        case 'exponential'
            Fk(:,w) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam);
            Fk(w,:) = exp([max(k,k(w(1))) + epsilon, max(k,k(w(2))) + epsilon]*gam)';
        end

    case 'deg-prod'
        
        switch modelform
        case 'powerlaw'
            Fk(:,w) = ([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam);
            Fk(w,:) = (([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon].^gam)');
        case 'exponential'
            Fk(:,w) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam);
            Fk(w,:) = exp([k*k(w(1)) + epsilon, k*k(w(2)) + epsilon]*gam)';
        end     

    end

	TopoData.k = k;
	TopoData.Fk = Fk;

    case 'neighbors'

	K = TopoData.K;
	Fk = TopoData.Fk;
	A = TopoData.A;
	u = TopoData.u;
	v = TopoData.v;

    uu = u(r);
    vv = v(r);
    x = A(uu,:);
    y = A(:,vv);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    K(uu,y) = K(uu,y) + 1;
    K(y,uu) = K(y,uu) + 1;
    K(vv,x) = K(vv,x) + 1;
    K(x,vv) = K(x,vv) + 1;
    switch modelform
        case 'powerlaw'
            Fk(uu,y) = (K(uu,y).^gam);
            Fk(y,uu) = Fk(uu,y)';
            Fk(vv,x) = (K(vv,x).^gam);
            Fk(x,vv) = Fk(vv,x)';
        case 'exponential'
            Fk(uu,y) = exp(K(uu,y)*gam);
            Fk(y,uu) = Fk(uu,y)';
            Fk(vv,x) = exp(K(vv,x)*gam);
            Fk(x,vv) = Fk(vv,x)';
    end

    TopoData.K = K;
	TopoData.Fk = Fk;
    TopoData.A = A;

    case 'matching'

	nei = TopoData.nei;
	degmat_sum = TopoData.degmat_sum;
	Fk = TopoData.Fk;
	A = TopoData.A;
	u = TopoData.u;
	v = TopoData.v;

	uu = u(r);
    vv = v(r);

    uu_nei = A(uu,:);
    vv_nei = A(vv,:);
    
    A(uu,vv) = 1;
    A(vv,uu) = 1;
        
    nei(uu,vv_nei) = nei(uu,vv_nei) + 1;
    nei(vv_nei,uu) = nei(vv_nei,uu) + 1;
    nei(vv,uu_nei) = nei(vv,uu_nei) + 1;
    nei(uu_nei,vv) = nei(uu_nei,vv) + 1;
        
    degmat_sum(uu,:) = degmat_sum(uu,:)+1;
    degmat_sum(vv,:) = degmat_sum(vv,:)+1;
    degmat_sum(:,uu) = degmat_sum(:,uu)+1;
    degmat_sum(:,vv) = degmat_sum(:,vv)+1;

    all_nei = [uu vv find(uu_nei) find(vv_nei)];
  
    switch modelform
        case 'powerlaw'
            %K_update = ((2 * nei(all_nei,:) ./ ( (degmat_sum(all_nei,:)<=2 & nei(all_nei,:)~=1)+(degmat_sum(all_nei,:) - (A(all_nei,:) * 2)) ) ) + epsilon);
            Fk_update = ( (2 * nei(all_nei,:) ./ ( (degmat_sum(all_nei,:)<=2 & nei(all_nei,:)~=1)+(degmat_sum(all_nei,:) - (A(all_nei,:) * 2)) ) ) + epsilon).^gam;
        case 'exponential'
            Fk_update = exp(( (2 * nei(all_nei,:) ./ ( (degmat_sum(all_nei,:)<=2 & nei(all_nei,:)~=1)+(degmat_sum(all_nei,:) - (A(all_nei,:) * 2)) ) ) + epsilon)*gam);
    end
    
    Fk(all_nei,:) = Fk_update;    
    Fk(:,all_nei) = Fk_update'; 

    TopoData.nei = nei;
	TopoData.degmat_sum = degmat_sum;
	TopoData.Fk = Fk;
    TopoData.A = A;

end

function TopoData = GetKSeed(A,modelform)

A = A>0;
n = length(A);
TopoData.A = A;
TopoData.k = sum(A,2);
k = TopoData.k;
c = clustering_coef_bu(A);

switch modelform
    
case 'clu-avg'        
        TopoData.c = clustering_coef_bu(A);
        TopoData.K = bsxfun(@plus,c(:,ones(1,n)),c')/2;
        
    case 'clu-diff'
        TopoData.c = clustering_coef_bu(A);
        TopoData.K = abs(bsxfun(@minus,c(:,ones(1,n)),c'));
        
    case 'clu-max'
        TopoData.c = clustering_coef_bu(A);
        TopoData.K = bsxfun(@max,c(:,ones(1,n)),c');
        
    case 'clu-min'
        TopoData.c = clustering_coef_bu(A);
        TopoData.K = bsxfun(@min,c(:,ones(1,n)),c');
        
    case 'clu-prod'
        TopoData.c = clustering_coef_bu(A);
        TopoData.K = c*c';
        
    case 'deg-avg'
        
        TopoData.K = bsxfun(@plus,k(:,ones(1,n)),k')/2;
        
    case 'deg-diff'

        TopoData.K = abs(bsxfun(@minus,k(:,ones(1,n)),k'));
        
    case 'deg-max'
        TopoData.K = bsxfun(@max,k(:,ones(1,n)),k');
        
    case 'deg-min'
        TopoData.K = bsxfun(@min,k(:,ones(1,n)),k');
        
    case 'deg-prod'
        TopoData.K = (k*k').*~eye(n);
        
    case 'neighbors'
        TopoData.K = (A*A).*~eye(n);
        
    case 'matching'
        TopoData.K = matching(A);
        TopoData.nei = (A*A).*~eye(n);
        degmat = repmat(TopoData.k',n,1);
        degmat_ = degmat';
        TopoData.degmat_sum = degmat + degmat_;
end