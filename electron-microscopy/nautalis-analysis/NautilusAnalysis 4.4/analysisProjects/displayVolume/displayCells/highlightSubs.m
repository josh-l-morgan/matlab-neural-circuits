function[blockSubs] = highlightSubs(subCell,blockSize);


%%
ys = blockSize(1);
xs = blockSize(2);
zs = blockSize(3);


obSub = subCell;
meanSub  = mean(obSub,1);


[y x z] = ind2sub([ys xs zs],1:(ys*xs*zs)); 

blockSubs = [y' x' z'];

blockSubs = blockSubs - repmat(mean(blockSubs,1),[size(blockSubs,1) 1]);
blockSubs = blockSubs + repmat(meanSub,[size(blockSubs,1) 1]);

