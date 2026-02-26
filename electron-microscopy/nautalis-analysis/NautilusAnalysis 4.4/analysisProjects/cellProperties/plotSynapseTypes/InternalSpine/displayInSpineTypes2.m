clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])

allEdges = obI.nameProps.edges(:,[2 1]);
seedList = [108 201 109 903 907];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);


inSpine = obI.nameProps.inSpine;
preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];
%inNum = [inSpine.tip];


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
clear histIn

for a = 1:length(preNames);
    
    getProps = (preID == preNames(a)) & (usePost>0);
    if sum(getProps)
        
        inList = [inSpine(getProps).in];
        inRat(a,:) = [sum(inList)/length(inList) sum(inList) length(inList)];
        ax.axID(a) = preNames(a);
        ax.synNum(a) = length(inList);
        ax.inList{a} = inList;
        ax.sumSpine(a) = sum(inList);
        ax.multNum(a) = sum(inList>1);
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

%%
clf
scatter(inRat(:,3),inRat(:,2))

for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    subplot(length(seedList),1,s)
    scatter(inRat(idx,3),inRat(idx,2));
    ylim([0 30])
    xlim([0 60])
    
end

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

%% Plot seed perspective
clf

for s = 1:length(seedList)

    preSeed = intersect(preNames,conTo(s).rgcList);
    linkedTC = intersect(postNames,conTo(s).tcrList);
    clear tempID sumSpine synNum multNum
    for t = 1:length(linkedTC)
        inList = [];
        for a = 1:length(preSeed)
            getProps = (postID == linkedTC(t)) & ...
                (preID == preSeed(a));
            if sum(getProps)
                inList = cat(2,inList,[inSpine(getProps).in]);
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

%% hist rat

minSyn = 5;
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
    percentGiantAx(s) = sum((useRat(:,1)>0.2))/size(useRat,1)*100
    
end

%% hist rat

minSyn = 10;
ratHistBin = [0:.01:.3];
clear percentGiant
for s = 1:length(seedList)
    [preSeed  idx]= intersect(preNames,conTo(s).preList);
    preRat = ax.multRat(idx);
    goodSyn = ax.synNum(idx)>=minSyn;
    useRat = preRat(goodSyn);
    histRat = histc(useRat,ratHistBin);
    subplot(length(seedList),1,s)
    bar(ratHistBin,histRat)
    %ylim([0 30])
    %xlim([0 1])
    pcsomething(s) = sum((useRat(:,1)>0.2))/size(useRat,1)*100
    
end

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
seedList = [108];
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
