

realBridgeDist = sort(sum(closeMat>0,1),'descend');


%% select links to use
useFrac = 1;
%%Put thumb on scale
realBridgeDist = realBridgeDist(1:round(length(realBridgeDist)*useFrac));

%realBridgeDist = realBridgeDist(2:end);
%realBridgeDist = realBridgeDist(realBridgeDist>0);

%% run simulation

totBridge = sum(realBridgeDist);
realCluster = rmsClust(realBridgeDist);
realLinked = sum(realBridgeDist>0)
cellNum = length(realBridgeDist);

reps = 10000;
randBridgeDists = zeros(reps,cellNum);
randCluster = zeros(reps,1);
randLinked = zeros(reps,1)
randBridgeDists(reps,cellNum);
for r = 1:reps
   
    randBridges = ceil(rand(totBridge,1)*cellNum);
    countBridges = hist(randBridges,1:cellNum);
    randCluster(r) =  cvrmse(countBridges);
    randLinked(r) =  sum(countBridges>0);

    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    
end

%% 
subplot(1,1,1)
binHist = [0:0.1:6]
binHistLinked = [0:1:cellNum];
histRandCluster = hist(randCluster,binHist);
histRandLinked = hist(randLinked,binHistLinked);

subplot(3,1,1)
bar(binHist,histRandCluster/max(histRandCluster(:)),'BarWidth',1)
hold on
scatter(realCluster,0,'r','LineWidth',5)
xlabel('Cluster size')
ylabel('Montecarlo distribution')
hold off
xlim([0 max(realCluster+realCluster/5,max(randCluster))])


subplot(3,1,2)

bar(binHistLinked,histRandLinked/max(histRandLinked(:)),'BarWidth',1)
hold on
scatter(realLinked,0,'r','LineWidth',5)
xlabel('LinkedNumber')
ylabel('Montecarlo distribution')
hold off
xlim([0 length(realBridgeDist)])


sortRandCluster = sort(randCluster,'ascend');
rand99 = [sortRandCluster(round(reps*.005)) sortRandCluster(end - round(reps*.005))] 
rand95 = [sortRandCluster(round(reps*.025)) sortRandCluster(end - round(reps*.025))] 

realCluster


sortRandLinked = sort(randLinked,'ascend');
rand99 = [sortRandLinked(round(reps*.005)) sortRandLinked(end - round(reps*.005))] 
rand95 = [sortRandLinked(round(reps*.025)) sortRandLinked(end - round(reps*.025))] 
meanRand = mean(sortRandLinked)
P = sum(randLinked<=realLinked)/reps

realLinked

sum(realBridgeDist)

%%
subplot(3,1,3)

meanRandBridgeDists = mean(randBridgeDists,1);

bar([ meanRandBridgeDists' realBridgeDist'])
xlabel('Cell bodies sorted by bridge number')
ylabel('Bridge number')
xlim([0 length(realBridgeDist)])


