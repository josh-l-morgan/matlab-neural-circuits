clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'butSize4.mat'])

allEdges = obI.nameProps.edges(:,[2 1]);
seedList = [108 201 903 907 109];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

%%Get internal spine data
inSpine = obI.nameProps.inSpine;
preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];
%inNum = [inSpine.tip];
hist([inSpine([inSpine.post] == 907).in])

allPres = unique(preID(preID>0));
allPres = intersect(allPres,useList.preList);
postNames = unique(postID(postID>0)); %unique(postID);
postNames = intersect(postNames,useList.postList);

preNames = [];
preCount = 0;
preSynNum = [];
for p = 1:length(allPres)
    %synNum = sum(allEdges(:,1) == allPres(p));
    synNum = sum(preID == allPres(p));
    if synNum>=0
        preCount = preCount + 1;
        preNames(preCount) = allPres(p);
        preSynNum(preCount) = synNum;
    end
end

%% Match anchors
voxLength = butSize.voxLength; % 0.15; %size of the voxels to be generated for measuring volumes
%anchorScale = [.0184 0.016 0.030]/voxLength;
anchorScale =  (  [obI.em.res].*[4 4 1])./ (1000*voxLength )  


obAnchors = obI.colStruc.anchors; % get object anchors
spineAnchors = obAnchors([inSpine.objectID],:);
spineAnchors = round(scaleSubs(double(spineAnchors),anchorScale));

spineEdges = [ [inSpine.post]' [inSpine.pre]'];

collectAnchors = [];

 


minDists = [];
inDiam = zeros(size(spineAnchors,1),1);
for a = 1:length(butSize.axList)
    
    
    butEdge = butSize.synDat(a).edges;
    
    butAnchors = butSize.synDat(a).anchors;
    goodAnchors = butAnchors>0;
    butAnchors(~goodAnchors) = 1;
    %butAnchors = round(scaleSubs(double(butAnchors),1./anchorScale));
    
    butVols = butSize.butVols{a};
    butVols = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2;
    
    collectAnchors = cat(1,collectAnchors,butAnchors);
    for s = 1:size(butAnchors,1)
        dists = sqrt((spineAnchors(:,1)-butAnchors(s,1)).^2 + ...
            (spineAnchors(:,2)-butAnchors(s,2)).^2 + (spineAnchors(:,3)-butAnchors(s,3)).^2);
        edgeFit = (spineEdges(:,1) == butEdge(s,1)) & ...
            (spineEdges(:,2) == butEdge(s,2)) ;
        
        minDistAll = min(dists);
        dists(~edgeFit) = 1000000;
                minDistAll = min(dists);
        minDist = min(dists);
        minDists = [minDists minDist];
        targ = find(dists == minDist,1);
        if ~isempty(targ) & goodAnchors(s) & (minDist < 5)
            
            inDiam(targ) = butVols(s);
            
        end
    end
end

subplot(2,1,1)
binDists = [0:1:10];
histDists = hist(minDists,binDists);
bar(binDists,histDists)

sum(inDiam>0)/length(inDiam)
subplot(2,1,2)
scatter(inDiam,inNum,'k')
scatter(inDiam,inNum,'k','filled')

%% scatter diam inNum by seed

clear percentGiant
for s = 1:length(seedList)
    
    [preSeed  idx] = intersect(preNames,conTo(s).rgcList);
    catDiam = [];
    catIn = [];
    hasPerf = preSeed*0;
    for a = 1:length(preSeed)
        getProp = find(preID == preSeed(a));
        useDiam = inDiam(getProp);
        useInNum = inNum(getProp);
        useInNum = useInNum(useDiam>0);
        useDiam = useDiam(useDiam>0);
    
        catDiam = [catDiam; useDiam];
        catIn = [catIn useInNum];
        hasPerf(a) = sum(useInNum);
    end
    
    pcPerf(s) = mean(hasPerf>0)*100;
    inFrac(s) = sum(catIn)/length(catIn)*100;
    pcInBut(s) = sum(catIn>0)/length(catIn)*100;
    subplot(length(seedList),1,s)
    scatter(catDiam,catIn,'k')
    ylim([0 11])
    xlim([0 4])

end
pcPerf
inFrac
pcInBut
%% scatter diam inNum by seed (seed only)

for s = 1:length(seedList)
    
     getProp = find(postID == seedList(s));
        useDiam = inDiam(getProp);
        useInNum = inNum(getProp);
        useInNum = useInNum(useDiam>0);
        useDiam = useDiam(useDiam>0);
        usePre = preID(getProp);
        useVol = diam2vol(useDiam);
        
    subplot(length(seedList),1,s)
    scatter(useDiam,useInNum,'k','filled')
       set(gca,'Ytick',[],'Xtick',[])

    ylim([0 11])
    xlim([0 4])

end



%%
inCount = [0:10];
butSizeBin = [0:.2:4]
for i = 1:length(inCount)
   subplot(length(inCount),1,i)
   histSize = hist(inDiam(inNum==inCount(end-i+1)),butSizeBin);
   bar(butSizeBin,histSize,.8,'k')
   ylim([0 200])
   set(gca,'Ytick',[],'Xtick',[])
    
    
end
     subplot(length(inCount),1,i)
   histSize = hist(inDiam(inNum==0),butSizeBin);
   bar(butSizeBin,histSize,.8,'k')
   ylim([0 200])
   set(gca,'Ytick',[],'Xtick',[])
    

%% Get cells

axList = preNames;
synapses = obI.nameProps.edges;
edges = synapses(:,1:2);

butSize.edges = edges;

rawSynAnchors = obI.colStruc.anchors(synapses(:,3),:);
synAnchors = scaleSubs(double(rawSynAnchors),anchorScale);




%% Filter post
usePost = postID * 0;
for i = 1:length(postID)
    
    %usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>200);
    usePost(i) =  sum(postNames == postID(i));
    usePre(i) = sum(preNames == preID(i));
end


%% Plot Pres
clf
plotNum = length(preNames);
boutBin = 0:10;
subY = floor(sqrt(length(preNames)));
subX = ceil(length(preNames)/subY);
inRat = zeros(length(preNames),3);

inBin = [0:12];
clear histIn ax

for a = 1:length(preNames);
    
    getProps = (preID == preNames(a)) & (usePost>0);
    if sum(getProps)
        
        inList = [inSpine(getProps).in];
        inRat(a,:) = [sum(inList)/length(inList) sum(inList) length(inList)];
        ax.axID(a) = preNames(a);
        ax.postIDs{a} = postID(getProps);
        ax.synNum(a) = length(inList);
        ax.inList{a} = inList;
        ax.sumSpine(a) = sum(inList);
        ax.multNum(a) = sum(inList>1);
        ax.anyNum(a) = sum(inList>0);
        ax.multRat(a) = sum(inList>1)/length(inList);
        ax.spineRat(a) = sum(inList)/length(inList);
        
        histIn(a,:) = hist(inList,inBin);
        
    end
    bar(inBin,sum(histIn,1))
    %subplot(plotNum,3,(a-1)*3+2);
    %     bar(boutBin,tipNum)
    %     xlim([min(tipNum)-1 max(tipNum)])
    %
    
    %    subplot(plotNum,3,(a-1)*3+3);
    %
    %     bar([filFrac 0])
    %     hold on
    %     bar([0 axAxFrac],'r')
    %     hold off
    %     ylim([0 1])
    
    
end
clf
scatter(inRat(:,3),inRat(:,2))

minSyn = 5;
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    
    getSynNum = ax.synNum(idx);
    getAnyNum = ax.anyNum(idx);
    
    pcAny = getAnyNum./getSynNum;
    pcAny = pcAny(getSynNum>=minSyn);
    meanAny(s) = mean(pcAny);
    stdAny(s) = std(pcAny);
    
    allPCany{s} = pcAny;
    subplot(length(seedList),1,s)
    scatter(inRat(idx,3),inRat(idx,2));
    ylim([0 30])
    xlim([0 60])
    
end

bigAny = [allPCany{[3 4 5]}];
meanBigAny = mean(bigAny);
stdBigAny = std(bigAny);
length(bigAny)

%% plot post

for a = 1:length(postNames);
    getProps = (postID == postNames(a)) & (usePre>0);
    if sum(getProps)
        
        inList = [inSpine(getProps).in];
        %inRat(a,:) = [sum(inList)/length(inList) sum(inList) length(inList)];
        tc.ID(a) = postNames(a);
        tc.synNum(a) = length(inList);
        tc.inList{a} = inList;
        tc.sumSpine(a) = sum(inList);
        tc.multNum(a) = sum(inList>1);
        tc.multRat(a) = sum(inList>1)/length(inList);
        tc.spineRat(a) = sum(inList)/length(inList);
        
        tcHistIn(a,:) = hist(inList,inBin);
        
    end
    bar(inBin,sum(tcHistIn,1))
    %subplot(plotNum,3,(a-1)*3+2);
    %     bar(boutBin,tipNum)
    %     xlim([min(tipNum)-1 max(tipNum)])
    %
    
    %    subplot(plotNum,3,(a-1)*3+3);
    %
    %     bar([filFrac 0])
    %     hold on
    %     bar([0 axAxFrac],'r')
    %     hold off
    %     ylim([0 1])
    
    
end

clf
for s = 1:length(seedList)
    [postSeed  idx]= intersect(postNames,conTo(s).tcrList);
    subplot(length(seedList)+1,1,s)
    scatter(tc.synNum(idx),tc.sumSpine(idx),'k');
    hold on
    targ = find(postNames == seedList(s));
    scatter(tc.synNum(targ),tc.sumSpine(targ),'r');
    
    ylim([0 70])
    %     xlim([0 60])
    
end



[postSeed  idx]= intersect(postNames,intersect(conTo(1).tcrList,conTo(2).tcrList));
subplot(length(seedList)+1,1,length(seedList)+1)
scatter(tc.synNum(idx),tc.sumSpine(idx),'k');
hold on
targ = find(postNames == seedList(s));
scatter(tc.synNum(targ),tc.sumSpine(targ),'r');

ylim([0 70])


%% Make getList mat

boutMat = [tc.ID(:) tc.spineRat(:) tc.synNum(:); ax.axID(:) ax.spineRat(:) ax.synNum(:)];

%% Make getList mat without seed cell
noSeedIn = zeros(length(ax.axID),1);
noSeedCount = zeros(length(ax.axID),1);
for a = 1:length(ax.axID)
   postA = ax.postIDs{a};
   inList = ax.inList{a};
   aSum = 0;
   aCount = 0;
   for p = 1:length(postA);
       if sum(seedList == postA(p)) == 0;
           aSum = aSum+inList(p);
           aCount = aCount+1;
       end
   end
   noSeedIn(a) = aSum;
   noSeedCount(a) = aCount;
    
end

boutMat = [tc.ID(:) tc.spineRat(:) tc.synNum(:); ...
    ax.axID(:) noSeedIn./noSeedCount noSeedCount];



%% Plot seed perspective
clf

for s = 1:length(seedList)

    preSeed = intersect(preNames,conTo(s).rgcList);
    linkedTC = intersect(postNames,conTo(s).tcrList);
    clear tempID sumSpine synNum multNum
    axCount = linkedTC*0;
    for t = 1:length(linkedTC)
        inList = [];
        for a = 1:length(preSeed)
            getProps = (postID == linkedTC(t)) & ...
                (preID == preSeed(a));
            if sum(getProps)
                inList = cat(2,inList,[inSpine(getProps).in]);
                            axCount(t) = axCount(t)+ 1;

            end
        end
        
        tempID(t) = linkedTC(t);
        sumSpine(t) = sum(inList);
        synNum(t) = length(inList);
        multNum(t) = sum(inList>1);
        %tcHistIn(a,:) = hist(inList,inBin);

    end
    
    subplot(length(seedList),1,s)
    scatter(synNum,sumSpine./synNum,20,'k','filled');
    hold on
    targ = find(tempID == seedList(s));
    scatter(synNum(targ),sumSpine(targ)./synNum(targ),40,'r','filled');
    plot([0 230],[sumSpine(targ)/synNum(targ) sumSpine(targ)/synNum(targ)],'r')
    ylim([0 3])
    xlim([0 230])
    

end


% 
% [postSeed  idx]= intersect(postNames,intersect(conTo(1).tcrList,conTo(2).tcrList));
% subplot(length(seedList)+1,1,length(seedList)+1)
% scatter(tc.synNum(idx),tc.sumSpine(idx),'k');
% hold on
% targ = find(postNames == seedList(s));
% scatter(tc.synNum(targ),tc.sumSpine(targ),'r');
% 
% ylim([0 70])



%% Plot seed vs hist of others
clf
minSyn = 10;
ratHistBin = [-.029:.03:.33];
for s = 1:length(seedList)

    preSeed = intersect(preNames,conTo(s).rgcList);
    linkedTC = intersect(postNames,conTo(s).tcrList);
    clear tempID sumSpine synNum multNum anyNum
    axCount = linkedTC*0;
    for t = 1:length(linkedTC)
        inList = [];
        for a = 1:length(preSeed)
            getProps = (postID == linkedTC(t)) & ...
                (preID == preSeed(a));
            if sum(getProps)
                inList = cat(2,inList,[inSpine(getProps).in]);
                            axCount(t) = axCount(t)+ 1;

            end
        end
        
        tempID(t) = linkedTC(t);
        sumSpine(t) = sum(inList);
        synNum(t) = length(inList);
        multNum(t) = sum(inList>1);
        anyNum(t) = sum(inList>0);

        
        %tcHistIn(a,:) = hist(inList,inBin);

    end
    
     useTC = find(synNum>=minSyn);

    subplot(length(seedList),2,s*2-1)
    %scatter(synNum,sumSpine./synNum,20,'k','filled');
     scatter(synNum,anyNum./synNum,20,'k','filled');

    hold on
    targ = find(tempID == seedList(s));
    useTC = setdiff(useTC,targ);
    %scatter(synNum(targ),sumSpine(targ)./synNum(targ),40,'r','filled');
    scatter(synNum(targ),anyNum(targ)./synNum(targ),40,'r','filled');

    %plot([0 230],[sumSpine(targ)/synNum(targ) sumSpine(targ)/synNum(targ)],'r')
    plot([0 230],[anyNum(targ)/synNum(targ) anyNum(targ)/synNum(targ)],'r')

    ylim([0 1])
    xlim([0 230])
    
        subplot(length(seedList),2,s*2)
        
    histTC = histc(anyNum(useTC)./synNum(useTC),ratHistBin)
    bar(ratHistBin,histTC,'k')
    hold on
        scatter(anyNum(targ)./synNum(targ),.01,40,'r','filled');
        xlim([-.1 .6])
    hold off
    
    

end


% 
% [postSeed  idx]= intersect(postNames,intersect(conTo(1).tcrList,conTo(2).tcrList));
% subplot(length(seedList)+1,1,length(seedList)+1)
% scatter(tc.synNum(idx),tc.sumSpine(idx),'k');
% hold on
% targ = find(postNames == seedList(s));
% scatter(tc.synNum(targ),tc.sumSpine(targ),'r');
% 
% ylim([0 70])


%% hist rat

minSyn = 10;
ratHistBin = [0:.1:1];
clear percentGiant
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    preRat = inRat(idx,:);
    goodSyn = preRat(:,3)>=minSyn;
    useRat = preRat(goodSyn,:);
    histRat = histc(useRat(:,1),ratHistBin);
    subplot(length(seedList),1,s)
    bar(ratHistBin,histRat,'k')
    %ylim([0 30])
    %xlim([0 1])
    	percentGiantAx(s) = sum((useRat(:,1)>0.1))/size(useRat,1)*100;
        [sum(useRat(:,1)> 0.1) size(useRat,1)]
        %percentGiantAx(s) = sum((useRat(:,2)>1))/size(useRat,1)*100

end 








%% compare on vs off seed for each axon
minSyn = 1;
ratHistBin = [0:.1:1];
clear percentGiant
clf
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    preRat = inRat(idx,:);
    
    preHist = preTo(allEdges,seedList(s));
    clear seedCon
    for p = 1:length(preSeed)
       seedCon(p,1) = preHist(preHist(:,1)==preSeed(p),2); 
    end
    
    goodSyn = preRat(:,3)>=minSyn;
    useRat = preRat(goodSyn,:);
    histRat = histc(useRat(:,1),ratHistBin);
    subplot(length(seedList),1,s)
        scatter(seedCon',preRat(:,2)./(preRat(:,3)-seedCon))

    %bar(ratHistBin,histRat,'k')
    %ylim([0 30])
    %xlim([0 1])
    	percentGiantAx(s) = sum((useRat(:,1)>0.1))/size(useRat,1)*100;
        [sum(useRat(:,1)> 0.1) size(useRat,1)]
        %percentGiantAx(s) = sum((useRat(:,2)>1))/size(useRat,1)*100

end
%% hist rat
%% ??
minSyn = 10;
ratHistBin = [-.029:.03:.33];
%ratHistBin = [0:.03:.3];

clear percentGiant
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    preRat = ax.multRat(idx);
    anyRat = ax.anyNum(idx)./ax.synNum(idx);
    goodSyn = ax.synNum(idx)>=minSyn;
    
    useRat = anyRat(goodSyn);
    histRat = histc(useRat,ratHistBin);
    subplot(length(seedList),1,s)
    bar(ratHistBin,histRat)
    %ylim([0 30])
    %xlim([0 1])
    pcsomething(s) = sum((useRat(:,1)>0.2))/size(useRat,1)*100
    
end

%% compare on vs off seed for each axon




%% seed bar
clf
minSyn = 10;
ratHistBin = [0:.1:1];
clear percentGiant seedHist
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    preRat = inRat(idx,:);
    goodSyn = preRat(:,3)>=minSyn;
    seedHist(s,:) = sum(histIn(idx,:),1);
    useRat = preRat(goodSyn,:);
    histRat = histc(useRat(:,1),ratHistBin);
    subplot(length(seedList),1,s)
    bar(ratHistBin,histRat)
    %ylim([0 30])
    %xlim([0 1])
    percentGiantBut(s) = sum(seedHist(s,3:end))/sum(seedHist(s,:))*100;
    bar(inBin,seedHist(s,:))
end
percentGiantBut


%% 3D conGraph


inSpine = obI.nameProps.inSpine;
preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

minCon = 10;
%seedList = [108];
useList108 = obI2cellList_seedInput_RGC_TCR(obI,seedList);

preNames = intersect(preID,useList108.preList);
postNames = intersect(postID,useList108.postList);

clear meanInPre meanInPost
for i = 1:length(preNames)
    isPre = preID == preNames(i);
    if sum(isPre)<minCon
        meanInPre(i) = -100;
    else
        meanInPre(i) = mean(inNum(isPre));
    end
end

for i = 1:length(postNames)
    isPost = postID == postNames(i);
    if sum(isPost)<minCon
        meanInPost(i) = -100;
    else
        meanInPost(i)  = mean(inNum(isPost));
        
    end
end

[sortedPreMean idx] = sort(meanInPre,'descend');
preNames = preNames(idx);
preNames = preNames(sortedPreMean>=0);
[sortedPostMean idx] = sort(meanInPost,'descend');
postNames = postNames(idx);
postNames = postNames(sortedPostMean>=0);

inCheck = [0:max(inNum)];
conI = zeros(length(preNames),length(postNames),length(inCheck));
for i = 1:length(preNames)
    for p = 1:length(postNames)
        for m = 1:length(inCheck);
            foundC = ((preID == preNames(i)) & (postID == postNames(p))...
                & (inNum == inCheck(m)));
            conI(i,p,m) = conI(i,p,m) + sum(foundC);
            
        end
    end
end

for i = 1:length(inCheck)
    subplot(length(inCheck),1,i);
    image(conI(:,:,i)*20)
end


conRat = conI(:,:,2)./conI(:,:,1);
clf
image(conRat*100)

clear colCon
colCon(:,:,1) = conI(:,:,1)*30;
colCon(:,:,2) = sum(conI(:,:,3:end),3)*200;
colCon(:,:,3) = sum(conI(:,:,5:end),3)*100;

image(uint8(colCon*1))


%% Only Big
clf
preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

minCon = 1;
preNames = unique(preID);
histPre = hist(preID,preNames)
preNames = preNames(histPre>=minCon);
preNames = preNames(preNames>0);
postNames = unique(postID);
histPost = hist(postID,postNames);
postNames = postNames(histPost>=minCon);
postNames = postNames(postNames>0);



clear meanInPre meanInPost
for i = 1:length(preNames)
    isPre = preID == preNames(i);
    meanInPre(i) = sum(inNum(isPre)>=2);
end

[sortedPreMean idx] = sort(meanInPre,'descend');
preNames = preNames(idx);
preNames = preNames(sortedPreMean>=1);

postNames = [];
for i = 1:length(preNames)
    isPre = preID == preNames(i);
    postNames = [postNames postID(isPre)];
end

postNames = unique(postNames(postNames>0));
for i = 1:length(postNames)
    isPost = postID == postNames(i);
    meanInPost(i)  = sum(inNum(isPost));
end

[sortedPostMean idx] = sort(meanInPost,'descend');
postNames = postNames(idx);

inCheck = [0:max(inNum)];
conI = zeros(length(preNames),length(postNames),length(inCheck));
for i = 1:length(preNames)
    for p = 1:length(postNames)
        for m = 1:length(inCheck);
            foundC = ((preID == preNames(i)) & (postID == postNames(p)) & (inNum == inCheck(m)));
            conI(i,p,m) = conI(i,p,m) + sum(foundC);
            
        end
    end
end
%
% for i = 1:length(inCheck)
%     subplot(length(inCheck),1,i);
%     image(conI(:,:,i)*20)
% end

%
% conRat = conI(:,:,2)./conI(:,:,1);
% clf
% image(conRat*100)

colCon  = zeros(size(conI,1),size(conI,2),3);
colCon(:,:,1) = conI(:,:,1)*50;
colCon(:,:,2) = sum(conI(:,:,3:end),3)*200;
colCon(:,:,3) = sum(conI(:,:,5:end),3)*0;

subplot(2,1,1)
image(uint8(colCon*1))
subplot(2,1,2)
someSyn = sum(conI,3);
noIn = conI(:,:,1)+ rand(size(conI,1),size(conI,2))*.2;
someIn = sum(conI(:,:,2:end),3) + rand(size(conI,1),size(conI,2))*.2;
scatter(noIn(someSyn>0),someIn(someSyn>0))

allBut = sum(sum(conI,3),2);
nonInBut = sum(conI(:,:,1),2)
manyInBut = sum(sum(conI(:,:,2:end),3),2)

inRat = manyInBut./allBut
mean(inRat)
std(inRat)
length(inRat)

mean(1-inRat)
std(1-inRat)
length(1-inRat)


%% type list



preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

minCon = 1;
preNames = unique(preID);
histPre = hist(preID,preNames)
preNames = preNames(histPre>=minCon);
preNames = preNames(preNames>0);
postNames = unique(postID);
histPost = hist(postID,postNames);
postNames = postNames(histPost>=minCon);
postNames = postNames(postNames>0);



for i = 1:length(preNames)
    foundC = ((preID == preNames(i))  & (inNum >= 2));
    preCount(i) =sum(foundC);
end


for p = 1:length(postNames)
    foundC = ( (postID == postNames(p)) & (inNum >=2));
    postCount(p) =sum(foundC);
end

isBig = [preNames(preCount>0) postNames(postCount>0)];


intersect(isBig, conTo(1).rgcList)


%% compare presynaptic bouton size to post synaptic internal spine
for i = 1:length(butSize.butDiam)
   medDiam(i) = median(butSize.butDiam{i}); 
   axSynNum(i) = length(butSize.butDiam{i}); 
end

medDiam = medDiam(axSynNum>5);
useAx = butSize.axList(axSynNum>5);
[medDiam idx] = sort(medDiam);
useAx = useAx(idx);

useTC = tc.ID(tc.synNum>5);
useIn = tc.spineRat(tc.synNum>5);

[useIn idx] = sort(useIn);
useTC = useTC(idx);


con = zeros(length(useAx),length(useTC));
inCon = con;
for a = 1:size(con,1)
    for p = 1:size(con,2)
        con(a,p) = sum((allEdges(:,1) == useAx(a)) & (allEdges(:,2)==useTC(p)));
        getProps = (postID == useTC(p)) &  (preID == useAx(a));
        if sum(getProps)
            inCon(a,p) = sum([inSpine(getProps).in]);
        end
    end
end

clf
hold on
cmap = [0 0 0; bluered(100)];
showcon = round((inCon./con)*100);
showcon(showcon>100) = 100;

for a = 1:size(con,1)
    for p = 1:size(con,2)
        if (con(a,p)>0) & ~(showcon(a,p))
        sc = scatter(useIn(p),medDiam(a),'markerfacecolor',cmap(showcon(a,p)+1,:),'MarkerEdgeColor','w');
        %set(sc,'edge',0)
        pause(.1)
        end
    end
end

for a = 1:size(con,1)
    for p = 1:size(con,2)
        if (con(a,p)>0) & (showcon(a,p))
        sc = scatter(useIn(p),medDiam(a),'markerfacecolor',cmap(showcon(a,p)+1,:),'MarkerEdgeColor','w');
        %set(sc,'edge',0)
        pause(.1)
        end
    end
end
hold off




useTC(useIn == .5)
useAx((medDiam >1.143) & (medDiam < 1.145))

useTC((useIn > .1428) & (useIn < .1430))


%% compare ax size to ax internal spine 

for i = 1:length(butSize.butDiam)
   medDiam(i) = median(butSize.butDiam{i}); 
   axSynNum(i) = length(butSize.butDiam{i}); 
end

medDiam = medDiam(axSynNum>5);
useAx = butSize.axList(axSynNum>5);

useAx2 = ax.axID(ax.synNum>10);
useIn = ax.spineRat(ax.synNum>10);

bothAx = intersect(useAx, useAx2);

for a = 1:length(bothAx)
   crossA(a,:) = [medDiam(useAx==bothAx(a)) useIn(useAx2==bothAx(a))] ;
end

scatter(jitter(crossA(:,1),.1),crossA(:,2),'k','filled')

%% compare tc size to tc internal spine 




useTC = tc.ID(tc.synNum>5);
useIn = tc.spineRat(tc.synNum>5);

clear tcSize
tcList = [];
for i = 1:length(useTC)
    tcSize{i} = [];
end
for i = 1:length(butSize.butDiam)
    edges = butSize.synDat(i).edges;
    diams = butSize.butDiam{i};
    for e = 1:size(edges,1)
       targ = find(useTC == edges(e,1));
       if ~isempty(targ)
           tcSize{targ} = [tcSize{targ} diams(e)];
       end
    end
    
end
clear tcMed tcLength
for t = 1:length(tcSize)
   tcMed(t) = median(tcSize{t}); 
   tcLength(t) = length(tcSize{t}); 
end

scatter(tcMed,useIn,'k','filled')
