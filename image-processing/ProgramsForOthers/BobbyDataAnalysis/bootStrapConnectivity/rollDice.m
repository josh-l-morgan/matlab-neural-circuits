
dNum = 2;
sideNum = 6;
reps = 10000;

res = fix(rand(reps,dNum)*sideNum)+1;
uniqueRes = unique(res,'rows');
numComb = size(uniqueRes,1)

histU = zeros(size(uniqueRes,1),sideNum);
for i = 1:size(uniqueRes,1)
   histU(i,:) = hist(uniqueRes(i,:),[1:1:sideNum]);
end


freqNum = hist(histU(:),[0:1:dNum]);
synNum = freqNum/sum(histU(:));
axNum = freqNum/numComb;
pairNum = freqNum/length(histU(:));

subplot(3,1,1)
bar([0:dNum],pairNum)

subplot(3,1,2)
bar([0:dNum],synNum)

subplot(3,1,3)
bar([0:dNum],axNum)

axNum(1:end-1)./axNum(2:end)

numAtMax = sum(histU(:)==dNum)/numComb






