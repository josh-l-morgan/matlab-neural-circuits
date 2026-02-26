pred2path(pred dist2seed)

vNum = length(pred);
uPred = unique(pred);
histPred = hist(pred,1:vNum);
onPath = pred>0;

noPred = onPath & ~histPred';
