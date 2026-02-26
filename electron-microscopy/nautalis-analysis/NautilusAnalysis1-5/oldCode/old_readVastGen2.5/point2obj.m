


TPN = [SPN(1:end-1) '_mat\'];
load([TPN 'sparseTiles.mat']);

allVals = cat(1,sparseTiles.id);
uniqueVals = unique(allVals);
countVals = hist(allVals,uniqueVals);
