function [Q,Ci]=modularity_Q(A,gamma)
%MODULARITY_UND     Optimal community structure and modularity
%
%   Ci = modularity_und(W);
%   [Ci Q] = modularity_und(W,gamma);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges.
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups.
%
%   Inputs:
%       W,
%           undirected weighted/binary connection matrix
%       gamma,
%           resolution parameter (optional)
%               gamma>1,        detects smaller modules
%               0<=gamma<1,     detects larger modules
%               gamma=1,        classic modularity (default)
%
%   Outputs:    
%       Ci,     optimal community structure
%       Q,      maximized modularity
%
%   Note:
%       This algorithm is essentially deterministic. The only potential
%       source of stochasticity occurs at the iterative finetuning step, in
%       the presence of non-unique optimal swaps. However, the present
%       implementation always makes the first available optimal swap and
%       is therefore deterministic.
%
%   References: 
%       Newman (2006) -- Phys Rev E 74:036104, PNAS 23:8577-8582.
%       Reichardt and Bornholdt (2006) Phys Rev E 74:016110.
%
%   2008-2016
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Dani Bassett, UCSB
%   Xindi Wang, Beijing Normal University
%   Roan LaPlante, Martinos Center for Biomedical Imaging

%   Modification History:
%   Jul 2008: Original (Mika Rubinov)
%   Oct 2008: Positive eigenvalues made insufficient for division (Jonathan Power)
%   Dec 2008: Fine-tuning made consistent with Newman's description (Jonathan Power)
%   Dec 2008: Fine-tuning vectorized (Mika Rubinov)
%   Sep 2010: Node identities permuted (Dani Bassett)
%   Dec 2013: Gamma resolution parameter included (Mika Rubinov)
%   Dec 2013: Detection of maximum real part of eigenvalues enforced (Mika Rubinov)
%               Thanks to Mason Porter and Jack Setford, University of Oxford
%   Dec 2015: Single moves during fine-tuning enforced (Xindi Wang)
%   Jan 2017: Removed node permutation and updated documentation (Roan LaPlante)

if ~exist('gamma','var')
    gamma = 1;
end

N=length(A);                            %number of vertices
% n_perm = randperm(N);                   %DB: randomly permute order of nodes
% A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Bg=B;
Ng=N;

while U(1)                              %examine community U(1)
    [V,D]=eig(Bg);
    [~,i1]=max(real(diag(D)));          %maximal positive (real part of) eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector

    S=ones(Ng,1);
    S(v1<0)=-1;
    q=S.'*Bg*S;                         %contribution to modularity

    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                         %maximal contribution to modularity
        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        while any(indg)                 %iterative fine-tuning
            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
            [qmax,imax]=max(Qit.*indg); %for i=1:Ng
            Sit(imax)=-Sit(imax);       %	Sit(i)=-Sit(i);
            indg(imax)=nan;             %	Qit(i)=Sit.'*Bg*Sit;
            if qmax>q                   %	Sit(i)=-Sit(i);
                q=qmax;                 %end
                S=Sit;
            end
        end

        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];                   %#ok<AGROW>
        end
    else                                %contribution nonpositive: U(1) is indivisible
        U(1)=[];
    end

    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    bg=B(ind,ind);
    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end

s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));
% Ci_corrected = zeros(N,1);              % DB: initialize Ci_corrected
% Ci_corrected(n_perm) = Ci;              % DB: return order of nodes to the order used at the input stage.
% Ci = Ci_corrected;                      % DB: output corrected community assignments
