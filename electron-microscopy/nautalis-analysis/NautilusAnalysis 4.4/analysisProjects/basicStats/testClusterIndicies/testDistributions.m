

N = 3;
reps = 1000;

zeroDist = zeros(1,N);
onesDist = ones(1,N);
biasedDist = zeroDist;
biasedDist(1) = 10;

[c r] = cvrmse(onesDist)
[c r] = cvrmse(biasedDist)



randNums = rand(reps,N*1);
randDist = zeros(reps,N+1);
for i = 1:reps
   randDist(i,:) = histc(randNums(i,:),[0:1/N:1]);
end
randDist = randDist(:,1:end-1);

clustIdx = zeros(reps,1);
for i = 1:size(randDist,1)
    clustIdx(i) = cvrmse(randDist(i,:));
end

subplot(2,1,1)
hist(randDist(:))
xlim([-.2 max(randDist)])

subplot(2,1,2)
hist(clustIdx)
xlim([-.2 max(clustIdx)])

sortRandCluster = sort(clustIdx,'ascend');
rand99 = [sortRandCluster(round(reps*.005)) sortRandCluster(end - round(reps*.005))] 
mean(clustIdx)

