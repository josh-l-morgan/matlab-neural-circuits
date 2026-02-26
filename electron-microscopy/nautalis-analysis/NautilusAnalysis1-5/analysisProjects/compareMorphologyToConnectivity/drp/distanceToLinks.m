
clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'cbDat.mat']);
load([MPN 'obI.mat']);
seedList = [108 201];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);

%% Get distances to all cell bodies from seed cells
cbID = cbDat.cbID;
cbCenters = cbDat.cbCenters;
cbDists = zeros(length(cbID),length(seedList));
for s = 1:length(seedList)
    targCB = find(cbID == seedList(s));

    cbDists(:,s) = sqrt((cbCenters(:,1)-cbCenters(targCB,1)).^2 + ...
    (cbCenters(:,2)-cbCenters(targCB,2)).^2 + ...
    (cbCenters(:,3)-cbCenters(targCB,3)).^2);

end

%% generate synapse matrix for all cb
synMat = useList.con;
cbSynMat = zeros(length(useList.preList),length(cbID));
for i = 1:length(cbID)
    if cbID(i)>0
        targ = find(useList.postList == cbID(i));
        if ~isempty(targ)
            cbSynMat(:,i) = synMat(:,targ);
        end
    end
end

%% nearest neighbor in axon
subplot(2,2,1)

allAD = [];
clear axDists
for a = 1:size(cbSynMat,1)
   axCon = cbSynMat(a,:); 
   postTo = find(axCon>0);
   postPos = cbCenters(postTo,:);
   c = 0;
   clear axDists
   for i = 1:length(postTo)
       for p = i+1:length(postTo)
           c = c+1;
                      
           aD(c) = sqrt((postPos(i,1)-postPos(p,1)).^2 + ...
               (postPos(i,2)-postPos(p,2)).^2+ ...
               (postPos(i,3)-postPos(p,3)).^2);
           
       end       
   end
   
   axDists{a} = aD;
   allAD = cat(1,allAD,aD');
    
end

histBinAx = [1:5:200];
histAxDists = hist(allAD,histBinAx);
bar(histBinAx,histAxDists)

%% nearest post
subplot(2,2,2)
postList = useList.postList;
usePost = intersect(postList,cbID);
usePost = setdiff(usePost,[108 109 201 903 907]);
con = useList.con;

postShareSyn = zeros(length(usePost));
postShareAx = postShareSyn;
postDistMat = zeros(length(usePost));
nearestCB2Post = zeros(length(usePost),1);
for i = 1:length(usePost)
    for p = 1:length(usePost)
        pos1 = cbCenters(find(cbID==usePost(i)),:);
        pos2 = cbCenters(find(cbID==usePost(p)),:);
        con1 = con(:,find(postList == usePost(i)));
        con2 = con(:,find(postList == usePost(p)));
        
        postShareAx(i,p) = sum(min(con1>0,con2>0));
        postShareSyn(i,p) = sum(min(con1,con2));
        postDistMat(i,p) = sqrt((pos1(1)-pos2(1)).^2 + ...
            (pos1(2)-pos2(2)).^2 + (pos1(3)-pos2(3)).^2);       
        
    end
    
    distToPost =   sqrt((pos1(1)-cbCenters(:,1)).^2 + ...
            (pos1(2)-cbCenters(:,2)).^2 + (pos1(3)-cbCenters(:,3)).^2);  
    distToPost = sort(distToPost,'ascend');
    nearestCB2Post(i) = distToPost(2);
end
        
shareRat = postShareSyn./postShareAx;
shareRat(isnan(shareRat)) = 0;
hist(shareRat(:))


plotComp = (postShareSyn>=0) & (postDistMat>0);
scatter(postDistMat(plotComp),postShareSyn(plotComp))

synBinWidth = 5;
histMeanSynBin = [0:synBinWidth:300];
meanSyn = zeros(length(histMeanSynBin)-1,1);
meanAx = meanSyn;
for i = 1:length(histMeanSynBin)-1
   inWin = (postDistMat>histMeanSynBin(i)) & (postDistMat<=histMeanSynBin(i+1));
   meanSyn(i) = mean(postShareSyn(inWin & (postShareSyn>0)));
   meanAx(i) = mean(postShareAx(inWin));
   meanRat(i) = mean(shareRat(inWin & (shareRat>0)));
end
hold on
plot(histMeanSynBin(2:end)-synBinWidth/2,meanSyn,'r','linewidth',4)
plot(histMeanSynBin(2:end)-synBinWidth/2,meanAx,'g','linewidth',3)
%plot(histMeanSynBin(2:end)-synBinWidth/2,meanSyn./meanAx,'m','linewidth',2)
plot(histMeanSynBin(2:end)-synBinWidth/2,meanRat,'k','linewidth',2)

hold off

%%
subplot(2,2,3)
meanSeedPos =(cbCenters(find(cbID == 108),:)+ cbCenters(find(cbID == 201),:))/2
distToMeanSeed =   sqrt((meanSeedPos(1)-cbCenters(:,1)).^2 + ...
            (meanSeedPos(2)-cbCenters(:,2)).^2 + (meanSeedPos(3)-cbCenters(:,3)).^2);  
[a idx] = intersect(cbID,usePost);
postDistToMeanSeed = distToMeanSeed(idx);
dist2CenterMat = repmat(postDistToMeanSeed,[1,length(postDistToMeanSeed)]);

meanDistThresh = 100;
histMinBin = [0:5:100];


maskDist = postDistMat;
maskDist(postShareSyn<1) = 10000;
maskDist(postDistMat == 0) = 10000;
%maskDist(dist2CenterMat> meanDistThresh) = 10000;
minDistConnected = min(maskDist,[],2);

histMinConnected = hist(minDistConnected(postDistToMeanSeed<meanDistThresh,:),histMinBin);
bar(histMinBin,histMinConnected)

maskDist = postDistMat;
maskDist(postShareSyn>0) = 10000;
maskDist(postDistMat == 0) = 10000;
%maskDist(dist2CenterMat> meanDistThresh) = 10000;
minDistDisconnected = min(maskDist,[],2);

histMinDisConnected = hist(minDistDisconnected(postDistToMeanSeed<meanDistThresh,:),histMinBin);



histMinCB = hist(nearestCB2Post(postDistToMeanSeed<meanDistThresh,:),histMinBin);

colormap([1 0 0; 0 1 0; 0 0 1])
histMinBoth = [histMinConnected' histMinDisConnected' histMinCB'];
bar(histMinBin,histMinBoth)
plot(histMinBin,histMinBoth)

hold off
plot(histMinBin,histMinCB,'b')
hold on
plot(histMinBin,histMinConnected,'g')
plot(histMinBin,histMinDisConnected,'r')
hold off

ranksum(minDistConnected,minDistDisconnected)

%% repeat using networks instead of single axons
subplot(2,2,4)

seedConMat1 = repmat(con(:,find(postList == seedList(1))),[1 size(con,2)])
seedConMat2 = repmat(con(:,find(postList == seedList(2))),[1 size(con,2)])
shareSeedCon1 = sum(min(seedConMat1,con),1);
shareSeedCon2 = sum(min(seedConMat2,con),1);
[a idx] = intersect(postList,usePost);
shareSeed1 = shareSeedCon1(idx);
shareSeed2 = shareSeedCon2(idx);

clear shareSeedMat1 shareSeedMat2
for i = 1:length(shareSeed1)
    shareSeedMat1(i,:) = shareSeed1 * shareSeed1(i);
    shareSeedMat2(i,:) = shareSeed2 * shareSeed2(i);

end
shareSeed = shareSeedMat1 + shareSeedMat2;




maskDist = postDistMat;
maskDist(shareSeed==0) = 10000;
maskDist(postDistMat == 0) = 10000;
%maskDist(dist2CenterMat> meanDistThresh) = 10000;
minDistConnected = min(maskDist,[],2);

histMinConnected = hist(minDistConnected(postDistToMeanSeed<meanDistThresh,:),histMinBin);
bar(histMinBin,histMinConnected)

maskDist = postDistMat;
maskDist(shareSeed>0) = 10000;
maskDist(postDistMat == 0) = 10000;
%maskDist(dist2CenterMat> meanDistThresh) = 10000;
minDistDisconnected = min(maskDist,[],2);

histMinDisConnected = hist(minDistDisconnected(postDistToMeanSeed<meanDistThresh,:),histMinBin);



histMinCB = hist(nearestCB2Post(postDistToMeanSeed<meanDistThresh,:),histMinBin);

colormap([1 0 0; 0 1 0; 0 0 1])
histMinBoth = [histMinConnected' histMinDisConnected' histMinCB'];
bar(histMinBin,histMinBoth)
plot(histMinBin,histMinBoth)

hold off
plot(histMinBin,histMinCB,'b')
hold on
plot(histMinBin,histMinConnected,'g')
plot(histMinBin,histMinDisConnected,'r')
hold off

ranksum(minDistConnected,minDistDisconnected)


%%
%{

allNearest = [];
clear axNearest
for a = 1:size(cbSynMat,2)
   axCon = cbSynMat(a,:); 
   postTo = find(axCon>0);
   postPos = cbCenters(postTo,:);
   c = 0;
   clear axDists
   for i = 1:length(postTo)
       for p = 1:length(postTo)
           
                      
           aD = sqrt((postPos(i,1)-postPos(p,1)).^2 + ...
               (postPos(i,2)-postPos(p,2)).^2+ ...
               (postPos(i,3)-postPos(p,3)).^2);
           matDists(i,p) = aD;
       end       
   end
   matDists(matDists == 0) = 10000;
   minDists = min(matDists,[],2)
      
   allNearest =  cat(1,allNearest,minDists);
   axNearest{a} = minDists;
    
end

histBinAx = [1:5:200];
histAxNearest = hist(allNearest,histBinAx);
bar(histBinAx,histAxNearest)



%% find links to seed
sharedSyn = zeros(length(cbID),length(seedList));
sharedAx = sharedSyn;
for s = 1:length(seedList)
    targCB = find(cbID == seedList(s));
    seedMat = repmat(cbSynMat(:,targCB),[1 length(cbID) ]);
    sharedSyn(:,s) = sum(min(cbSynMat,seedMat),1)';
    sharedAx(:,s) = sum(min(cbSynMat,seedMat)>0,1)'   ; 
end
isLinked = sharedAx>0;

%% Bin distances
binSize = 10;
binRad = binSize/2;
binRange = [binRad:binSize:100];
linkedAt = zeros(length(binRange),length(seedList));
totalAt = linkedAt;
clear groupCells
for i = 1:length(binRange);
    checkPos = binRange(i);
    for s = 1:length(seedList)
        inRange = (cbDists(:,s)>(checkPos-binRad))& (cbDists(:,s)<=(checkPos+binRad));
        linkedAt(i,s) = sum(isLinked(inRange,s));
        totalAt(i,s) = sum(inRange);
        groupCells{i,s}  = find(inRange);
    end
end
    
for s = 1:length(seedList)
    subplot(3,1,s)
    bar(binRange,totalAt(:,s),'b')
    hold on
   bar(binRange,linkedAt(:,s),'r')
   hold off
end

   subplot(3,1,3)
   plot(binRange,linkedAt./totalAt)

   %% Redistribute links
   
   recreateLinks = isLinked * 0;
   for b = 1:size(groupCells,1)
       for s = 1:length(seedList)
           if linkedAt(b,s)
               possible = groupCells{b,s};
               pick  = randperm(length(possible),linkedAt(b,s));
               picked = possible(pick);
               recreateLinks(possible,s) = isLinked(possible,s);
           end
       end
   end
   
   sum(recreateLinks);
   realSame = sum(sum(recreateLinks,2) == 2);
   
   
   reps = 1000;
   
   for r = 1:reps;
   randLinked = isLinked * 0;
   for b = 1:size(groupCells,1)
       for s = 1:length(seedList)
           if linkedAt(b,s)
               possible = groupCells{b,s};
               pick  = randperm(length(possible),linkedAt(b,s));
               picked = possible(pick);
               randLinked(picked,s) = 1;
           end
       end
   end
   
   sum(randLinked);
   randSame(r) = sum(sum(randLinked,2) == 2);
   end
   
   
  mcHistBin = [0:1:30];
  histSame = hist(randSame,mcHistBin);
  bar(mcHistBin,histSame,'y')
  hold on
  scatter(realSame,10,40,'r','filled')
  hold off
  
  
   
%}



