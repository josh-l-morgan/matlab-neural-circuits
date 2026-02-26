%%How many 0 bridge cells need to be thrown away in order for the
%%distribution of bridges to resemble the observed distribution

testFracs = [1];
trackMeanRandCluster = zeros(1,length(testFracs));
reps = 1000;
realBridgeDistRaw = sort(closeBridge,'descend');



for t = 1:length(testFracs)
    
realBridgeDist = sort(closeBridge,'descend');
useFrac = testFracs(t);
realBridgeDist = realBridgeDist(1:round(length(realBridgeDist)*useFrac));

totBridge = sum(realBridgeDist);
realCluster = rmsClust(realBridgeDist);
cellNum = length(realBridgeDist);

randBridgeDists = zeros(reps,cellNum);
randCluster = zeros(reps,1);
for r = 1:reps
   
    randBridges = ceil(rand(totBridge,1)*cellNum);
    uBridges = unique(randBridges);
    countBridges = hist(randBridges,uBridges);
    %randCluster(r) = sum(countBridges.^2)/totBridge;
    randCluster(r) =  rmsClust(countBridges);
    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    
end

%% 
subplot(1,1,1)
binHist = [0:0.05:6];
histRandCluster = hist(randCluster,binHist);
bar(binHist,histRandCluster/max(histRandCluster(:)),'BarWidth',1);
hold on
scatter(realCluster,0,'r','LineWidth',5)
xlabel('Cluster size')
ylabel('Montecarlo distribution')
hold off
xlim([0 max(realCluster+realCluster/5,max(randCluster))])

%%

meanRandBridgeDists = mean(randBridgeDists,1);

bar([ meanRandBridgeDists' realBridgeDist'])
xlabel('Cell bodies sorted by bridge number')
ylabel('Bridge number')

trackMeanRandCluster(t) = mean(randCluster);
trackRealCluster(t) = realCluster;
trackClusterDifference(t) = mean(randCluster)-realCluster;

barDif(t) = sum(abs(realBridgeDist-meanRandBridgeDists))/totBridge;
pause(1)
end

plot(testFracs,trackMeanRandCluster)
hold on
plot(testFracs,trackRealCluster)
plot(testFracs,abs(trackClusterDifference))
hold off

plot(testFracs,barDif);
plot(testFracs,(trackMeanRandCluster-trackRealCluster))

clusterDifference = abs(trackMeanRandCluster-realCluster);
%bestFit = find(





%%  Run again with defined fraction of real

realBridgeDist = sort(closeBridge,'descend');
realBridgeDist = realBridgeDist(1:round(length(closeBridge)/3));

totBridge = sum(realBridgeDist);
realCluster = rmsClust(realBridgeDist);
cellNum = length(realBridgeDist);

randBridgeDists = zeros(reps,cellNum);
randCluster = zeros(reps,1);
for r = 1:reps
   
    randBridges = ceil(rand(totBridge,1)*cellNum);
    uBridges = unique(randBridges);
    countBridges = hist(randBridges,uBridges);
    %randCluster(r) = sum(countBridges.^2)/totBridge;
    randCluster(r) =  rmsClust(countBridges);
    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    
end

%% 
subplot(1,1,1)
binHist = [0:0.05:6];
histRandCluster = hist(randCluster,binHist);
bar(binHist,histRandCluster/max(histRandCluster(:)),'BarWidth',1);
hold on
scatter(realCluster,0,'r','LineWidth',5)
xlabel('Cluster size')
ylabel('Montecarlo distribution')
hold off
xlim([0 max(realCluster+realCluster/5,max(randCluster))])

%%

meanRandBridgeDistsShort = mean(randBridgeDists,1);
meanRandBridgeDistsShortBuff = meanRandBridgeDists*0;
meanRandBridgeDistsShortBuff(1:length(meanRandBridgeDistsShort)) = meanRandBridgeDistsShort;


bar([ realBridgeDistRaw' meanRandBridgeDists'  meanRandBridgeDistsShortBuff'])
xlabel('Cell bodies sorted by bridge number')
ylabel('Bridge number')

trackMeanRandCluster(t) = mean(randCluster);
trackRealCluster(t) = realCluster;
trackClusterDifference(t) = mean(randCluster)-realCluster;

barDif(t) = sum(abs(realBridgeDist-meanRandBridgeDists))/totBridge;
pause(1)



