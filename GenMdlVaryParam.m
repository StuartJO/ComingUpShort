function [adj, b, Edgeprob] = GenMdlVaryParam(PDMs,params,m,Input) 

% load('scha200_iFOD2_thr.mat')
% load('PDMs_Scha200.mat', 'CGE_dist')
% PDMs{1} = A_dist; PDMs{2} = CGE_dist;
% params= [-.1 1 1 5];
% PDMfunc = {'exponential','powerlaw'};
% m = 650;
% normType = 'max';
% AddMult = 'Add';

PDMfunc = Input.PDMfunc;
Timing = Input.Timing;
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

eta = params(1,1);
if Timing == 1

TimeSteps = length((mseed + 1):m);

eta_vals = linspace(0,eta,TimeSteps);

PDMparams(1) = eta_vals(1);

elseif Timing == 2

TimeSteps = length((mseed + 1):m);

eta_vals = linspace(-eta,eta,TimeSteps);

PDMparams(1) = eta_vals(1);

end

a = ones(1,length(indx));

b = zeros(m,1);

nn = length(indx);

nPDMParams = length(PDMs);

d = zeros(nPDMParams,nn);
Df = zeros(nPDMParams,nn);
PDf = zeros(nPDMParams,nn);

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

Edgeprob = zeros(Nsteps,nn);

for i = 1:Nsteps

EdgeInd = EdgeInds(i);

Edgeprob(i,:) = P;

    C = [0 cumsum(P)];
    r = sum(rand*C(end) >= C);
    b(EdgeInd) = r;
    
    if r<= 0
       save('DEBUG.mat') 
    end
    a(r) = 0;
      if Timing > 0  
    switch PDMfunc{1}
        case 'powerlaw'
            Df(1,:) = d(1,:).^eta_vals(i);
        case 'exponential'       
            Df(1,:) = exp(eta_vals(i)*(d(1,:)));    
    end
      end
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
%b_ = b;
b = indx(b);

adj = zeros(n);
adj(b) = 1;
adj = adj + adj';
% 
% c = zeros(nn,1);
% 
% c(b_) = length(b):-1:1;
% 
% scatter(d(1,:),d(2,:),30,c,'filled')
% 
% for i = 1:m
%    scatter(d(1,:),d(2,:),30,log(prod(ProbATTime{i},1)),'filled') 
%    pause(.1) 
% end


function TopoData = CalcTopoForm(modeltype,modelform,r,TopoData,gam,epsilon)

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
  
        %save('DEBUG.mat')
    
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