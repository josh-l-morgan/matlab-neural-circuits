

totBridge = sum(closeBridge);
realCluster = sum(closeBridge.^2)/totBridge;
realBridgeDist = sort(closeBridge,'descend');
cellNum = length(closeBridge);

reps = 100000;
randBridgeDists = zeros(reps,cellNum);
randCluster = zeros(reps,1)
for r = 1:reps
   
    randBridges = ceil(rand(totBridge,1)*cellNum);
    uBridges = unique(randBridges);
    countBridges = hist(randBridges,uBridges);
    %randCluster(r) = sum(countBridges.^2)/totBridge;
    randCluster(r) =  rms((countBridges-mean(countBridges))/mean(countBridges));;
    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    
end

%% 
subplot(1,1,1)
binHist = [0:0.05:6]
histRandCluster = hist(randCluster,binHist);
bar(binHist,histRandCluster/max(histRandCluster(:)),'BarWidth',1)
hold on
scatter(realCluster,0,'r','LineWidth',5)
xlabel('Cluster size')
ylabel('Montecarlo distribution')
hold off
xlim([0 5])

%%

meanRandBridgeDists = mean(randBridgeDists,1);

bar([ meanRandBridgeDists' realBridgeDist'])
xlabel('Cell bodies sorted by bridge number')
ylabel('Bridge number')
