realMat = closeMat %from connectVScb
realMat = realMat(sum(realMat,2)>0,:);

realSyn = sum(realMat,1);
realLinks = sum(realMat>0,1);
sortLinks = sort(realLinks,'descend');
realSpread = sum(realLinks>0);


realCluster = rmsClust(sortLinks);
realLinked = sum(sortLinks>0);

sumRealSyn = sum(realMat,1);
realSynCluster = cvrmse(sumRealSyn);



%% Link vs mean Syn
wasLink = realLinks(realLinks>0);
wasSyn = realSyn(realLinks>0);
realMeanSyn = wasSyn./wasLink;
scatter(wasLink,realMeanSyn)
[sortWasLink idx] = sort(wasLink,'descend');
sortMeanSyn = realMeanSyn(idx);

%%
tic
reps = 100000;
axNum = size(realMat,1);
closeNum = size(realMat,2);
totBridge = sum(sortLinks);
randLinked = zeros(reps,1);
randBridgeDists = zeros(reps,closeNum);
randCluster = zeros(reps,1);
newUse = useList;
%realMod = use2Modularity(useList);
for r = 1:reps
    
    newMat = realMat * 0;
    for a = 1:axNum
        newMat(a,:) = realMat(a,randperm(closeNum));
    end
    randLinked(r) =sum(sum(newMat,1)>0); 
    %image(newMat*20),pause(.3)
    
    %newUse.con = newMat;
    %randMod(r) = use2Modularity(newUse);
    
    
    sumRandSyn =  sum(newMat,1);
    randBridges = ceil(rand(totBridge,1)*closeNum);
    countBridges = hist(randBridges,1:closeNum);
    randCluster(r) =  cvrmse(countBridges);
    randLinked(r) =  sum(countBridges>0);
    randSynCluster(r) = cvrmse(sumRandSyn);
    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    sumRandSyns(r,:) = sumRandSyn;
    
    
    
%     
%     
end


toc






%% Convert to previous
%sortLinks = sort(closeBridge,'descend');


%{
%% select links to use
useFrac = 1;
%%Put thumb on scale
sortLinks = sortLinks(1:round(length(sortLinks)*useFrac));

%sortLinks = sortLinks(2:end);
%sortLinks = sortLinks(sortLinks>0);

%% run simulation

totBridge = sum(sortLinks);
realCluster = rmsClust(sortLinks);
realLinked = sum(sortLinks>0)
closeNum = length(sortLinks);

reps = 1000000;
randBridgeDists = zeros(reps,closeNum);
randCluster = zeros(reps,1);
for r = 1:reps
   
    randBridges = ceil(rand(totBridge,1)*closeNum);
    countBridges = hist(randBridges,1:closeNum);
    randCluster(r) =  cvrmse(countBridges);
    randLinked(r) =  sum(countBridges>0);

    randBridgeDists(r,1:length(countBridges)) = sort(countBridges,'descend');
    
end
%}
%% 
subplot(1,1,1)
binHist = [0:0.1:6]
binHistLinked = [0:1:closeNum];
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
xlim([0 max(realLinked+realLinked/5,max(randLinked))])


sortRandCluster = sort(randCluster,'ascend');
rand99 = [sortRandCluster(round(reps*.005)) sortRandCluster(end - round(reps*.005))] 
rand95 = [sortRandCluster(round(reps*.025)) sortRandCluster(end - round(reps*.025))] 

realCluster


sortRandLinked = sort(randLinked,'ascend');
rand99 = [sortRandLinked(round(reps*.005)) sortRandLinked(end - round(reps*.005))] 
rand95 = [sortRandLinked(round(reps*.025)) sortRandLinked(end - round(reps*.025))] 
mean(sortRandLinked)
P = sum(randLinked<=realLinked)/reps

realLinked

%% randSynCluster
maxSynClust = max(max(randSynCluster),realSynCluster);
binHist2 = [0:maxSynClust/20:maxSynClust];
subplot(3,1,1)
histRandSynCluster = hist(randSynCluster,binHist2);
bar(binHist2,histRandSynCluster/max(histRandSynCluster(:)),'BarWidth',1)
hold on
scatter(realSynCluster,0,'r','LineWidth',5)
xlabel('Cluster size')
ylabel('Montecarlo distribution')
hold off
xlim([0 maxSynClust])



%%
subplot(3,1,3)

meanRandBridgeDists = mean(randBridgeDists,1);

bar([ meanRandBridgeDists' sortLinks'])
xlabel('Cell bodies sorted by bridge number')
ylabel('Bridge number')


