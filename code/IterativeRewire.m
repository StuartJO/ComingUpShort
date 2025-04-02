function [FitMeasures,EdgeOverlap,DistRewires] = IterativeRewire(A,D,varargin)

p = inputParser;

addRequired(p,'A',@(x) ismatrix(x))
addRequired(p,'D',@(x) ismatrix(x))

totalEdges = nnz(A)/2;

n = length(A);

defaultRepeats=100;

DistMatchType = 'none';
defaultOrder='random';

expectedOrders = {'random','shortest','longest','any'};
DistMatchTypes = {'none','exact','prob','invprob'};

deaultRepeatRewires = false;

addParameter(p,'Rewires',totalEdges,@(x) isscalar(x) & ((x>0 & x<1) | rem(x,1)==0))
addParameter(p,'Repeats',defaultRepeats,@(x) rem(x,1)==0 & isscalar(x))
addParameter(p,'DistMatchType',DistMatchType,@(x) any(validatestring(x,DistMatchTypes)))
addParameter(p,'RewireOrder',defaultOrder, @(x) any(validatestring(x,expectedOrders)))
addParameter(p,'RepeatRewires',deaultRepeatRewires,@(x) islogical(x) & isscalar(x))

parse(p,A,D,varargin{:});

NEdges = p.Results.Rewires;
DistMatchType = p.Results.DistMatchType;
repeats = p.Results.Repeats;
RewireOrder = p.Results.RewireOrder;
RepeatRewires = p.Results.RepeatRewires;

switch DistMatchType
    case {'exact','prob','invprob'}
DistMatch = true;
    otherwise
DistMatch = false;
end

if NEdges < 1 && NEdges > 0
    numIterations = round(totalEdges*NEdges);
else
    numIterations = NEdges;
end

if RepeatRewires
    A_ = zeros(n);
else
    A_ = A;
end

% Initialize the number of iterations

I = eye(n);

II = triu(ones(n));

A_vals{1} = sum(A,2);
A_vals{2} = clustering_coef_bu(A);
A_vals{3} = betweenness_bin(A)';
A_vals{4} = D(triu(A,1) > 0);

ANetStats = CalcNetStats4Eval(A,D,A_vals);

origEdges = find(triu(A,1));
origEdgesD = D(origEdges);
switch RewireOrder
    case 'shortest'
        [~,ord] = sort(origEdgesD,'ascend');
    case 'longest'
        [~,ord] = sort(origEdgesD,'descend');
end

switch RewireOrder
    case {'shortest','longest'}
        if repeats>1 && strcmp(DistMatchType,'exact')
            repeats = 1;
            warning('When rewiring according to (exact) distances and trying to distance matching, the algorithm is deterministic. Setting ''Repeats'' to 1')
        end
end

FitMeasures.maxKS = zeros(repeats,numIterations);
FitMeasures.TND = zeros(repeats,numIterations);
FitMeasures.maxRMSE = zeros(repeats,numIterations);
FitMeasures.TFdiff = zeros(repeats,numIterations);
FitMeasures.maxRd = zeros(repeats,numIterations);

EdgeOverlap = zeros(repeats,numIterations);
FitMeasures.KS = zeros(repeats,numIterations,4);

DistRewires = zeros(numIterations,2,repeats);

h = waitbar(0);
h2 = waitbar(0);
for repeat = 1:repeats
B = A;
switch RewireOrder
    case 'random' 
    ord = randperm(length(origEdges));
end
    for iteration = 1:numIterations        
        switch RewireOrder
            case {'shortest','longest','random'}
                existingEdgeIdx = origEdges(ord(iteration));  
                [row_existing, col_existing] = ind2sub(n,existingEdgeIdx);
            case 'any'
                existingEdges = find(B);
                ind = randi(length(existingEdges));
                existingEdgeIdx = existingEdges(ind);
                [row_existing, col_existing] = ind2sub(n,existingEdgeIdx);
                B(row_existing, col_existing) = 0;
                B(col_existing, row_existing) = 0;
        end
        
    if DistMatch

        existingEdge_edgeDist = D(existingEdgeIdx);

        % Generate a list of non-existing edges
        
        nonexistingIdx = find( (B+A_+II) == 0);

        nonexisting_dist = D(nonexistingIdx);

        EdgeDistdiff = abs(existingEdge_edgeDist-nonexisting_dist);

        switch DistMatchType
            case 'exact'
                [~,Toswap] = min(EdgeDistdiff);
            case 'prob'
                wei = EdgeDistdiff.^-.01;
                [~,Toswap] = datasample(nonexistingIdx,1,'Weights',wei./max(wei));
            case 'invprob'
                wei = EdgeDistdiff.^2;
                [~,Toswap] = datasample(nonexistingIdx,1,'Weights',wei./max(wei));
        end

        [row_nonexisting, col_nonexisting] = ind2sub(n,nonexistingIdx(Toswap));

    else
        nonexistingIdx = find( (B+A_+I) == 0);
        Toswap = randi(length(nonexistingIdx));
        [row_nonexisting, col_nonexisting] = ind2sub(n,nonexistingIdx(Toswap));
    end

        % Replace the existing edge with the non-existing edge
        B(row_existing, col_existing) = 0;
        B(col_existing, row_existing) = 0;
        B(row_nonexisting, col_nonexisting) = 1;
        B(col_nonexisting, row_nonexisting) = 1;

        DistRewires(iteration,:,repeat) = [D(row_existing, col_existing) D(row_nonexisting, col_nonexisting)];

        BNetStats = CalcNetStats4Eval(B,D);

        [FitMeasures.maxKS(repeat,iteration),~,FitMeasures.TND(repeat,iteration),FitMeasures.maxRMSE(repeat,iteration),FitMeasures.maxRd(repeat,iteration),FitMeasures.TFdiff(repeat,iteration),~,TopoCorrs,~,FitMeasures.KS(repeat,iteration,:)] = CalcEvalStats(A,B,ANetStats,BNetStats,0);
        FitMeasures.DegCorr(repeat,iteration) = TopoCorrs(1);
        AB = B.*A;
        EdgeOverlap(repeat,iteration) = (nnz(AB)/2)/totalEdges;

        waitbar(iteration / numIterations, h2, ['Finished ', num2str(iteration), '/', num2str(numIterations)]);
    end

    waitbar(repeat / repeats, h, ['Finished ', num2str(repeat), '/', num2str(repeats)]);
end