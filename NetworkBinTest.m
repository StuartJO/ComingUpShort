function Outcome = NetworkBinTest(A,B)

Avec = triu2vec(A)>0;
Bvec = triu2vec(B)>0;

TP = sum(Avec&Bvec);

TN = sum((Avec==0)&(Bvec==0));

FP = sum(Avec==0 & Bvec==1);

FN = sum(Avec==1 & Bvec==0);

FPR = FP/(FP+TN);

FNR = FN/(FN+TP);

SEN = TP/(TP+FN);

SPEC = TN/(TN+TP);

ACC = (TP+TN)/(TP+TN+FP+FN);

Outcome = [FPR FNR SEN SPEC ACC];
