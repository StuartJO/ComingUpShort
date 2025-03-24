function NetStats4Eval = CalcNetStats4Eval(A,A_dist,Avals)

if nargin < 3
Clu = clustering_coef_bu(A);
Bet = betweenness_bin(A)';
Deg = sum(A,2);
else
Deg = Avals{1};
Clu = Avals{2};
Bet = Avals{3};
end

NetStats4Eval.EdgeDists = A_dist(triu(A,1) > 0);

NodeDist = sum(A_dist.*A,2)./Deg;
NodeDist(isnan(NodeDist))=0;

Clo = closeness_bin(A,1)';

matchingA = mean(matching(A),2);

NetStats4Eval.NodeMeasures = [Deg Clu Clo NodeDist Bet matchingA];

NetStats4Eval.TF = corr(NetStats4Eval.NodeMeasures);

feature_sum = sum(NetStats4Eval.NodeMeasures);

NetStats4Eval.TF_ismissing = feature_sum==0;

Comps = graphComponents(A);
LargestComp = mode(Comps);
NetStats4Eval.GlobalTopo = [efficiency_bin(A) diffusion_efficiency(A(Comps==LargestComp,Comps==LargestComp)) modularity_Q(A) transitivity_bu(A)];

NetStats4Eval.meanNodeMeasures = mean(NetStats4Eval.NodeMeasures);
% NetStats4Eval.Clu = Clu;
% NetStats4Eval.Bet = Bet;
% NetStats4Eval.NodeDist = NodeDist;
% NetStats4Eval.Clo = Clo;