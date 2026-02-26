clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])

allEdges = obI.nameProps.edges(:,[2 1]);
targCells = [108 201 109 903 907];
conTo = makeConTo(obI,targCells);
isLinked = cat(2,conTo([3 4 5]).tcrList)
isLinked = unique(isLinked((isLinked>0) & (isLinked<1000)));

inSpine = obI.nameProps.inSpine;
preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

allPres = unique(preID(preID>0));

postNames = unique(postID(postID>0)); %unique(postID);
%postNames = unique(postID);

tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
tcrList = [conTo.tcrList];
postNames = intersect(postNames,tcrList);

preNames = [];
preCount = 0;
for p = 1:length(allPres)
    synNum = sum(allEdges(:,1) == allPres(p));
    if synNum>=10
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
end


%% Plot Pres
clf
plotNum = length(preNames)
boutBin = 0:10;
subY = floor(sqrt(length(preNames)));
subX = ceil(length(preNames)/subY);
for a = 1:plotNum;
    
    getProps = (preID == preNames(a)) & (usePost>0);
    inNum = hist([inSpine(getProps).in],boutBin);
    tipNum = hist([inSpine(getProps).tip],boutBin);
    showNum = tipNum
    showNum = showNum/preSynNum(a);
    
    subplot(subY,subX,a)
    bar(boutBin(1:end),showNum(1:end))
    %text(5,.5,num2str(preNames(a)))
    title(num2str(preNames(a)))
    
    %subplot(plotNum,3,(a-1)*3+1);
    xlim([min(showNum)-1 max(boutBin)])
    ylim([0 1])
    
    
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





% 
% 
% %% Plot Pres
% plotNum = length(preNames)
% boutBin = 0:10;
% for a = 1:plotNum;
%     
%     getProps = (preID == preNames(a)) & usePost;
%     getAllPre = (preID == preNames(a));
%     linNum = sum([synProp(getAllPre).lin]);
%     tcrNum = sum([synProp(getAllPre).tcr]);
%     inhibFrac(a) = linNum/(linNum+tcrNum);
%     
%     filFrac(a) = mean([synProp(getProps).filimentous]);
%     sumSyn(a) = sum(getProps);
%     %axAxFrac = mean([synProp(getProps).axoAxonic]);
%     
%     meanNear(a) = mean([synProp(getProps).boutonNeighbor]);
%     stdNear(a) = std([synProp(getProps).boutonNeighbor]);%/sqrt(sum(getProps));
%     meanSpin(a) = mean([synProp(getProps).spineNum]);
%     stdSpin(a) = std([synProp(getProps).spineNum]);%/sqrt(sum(getProps));
%     
%     multSpin(a) = mean([synProp(getProps).spineNum]>1);
%     multNear(a) =  mean([synProp(getProps).boutonNeighbor]>1);
%     
% end
% 
% 
% 
% subplot(plotNum,3,3);
% bar(filFrac)
% 
% % title('filFrac')
% 
% 
% subplot(plotNum,3,6);
% bar(inhibFrac)
% 
% % title('inhibFrac')
% 
% subplot(plotNum,3,9);
% bar(sumSyn)
% 
%  title('sumSyn')
% 
% subplot(plotNum,3,12);
% bar(meanNear)
% 
%  title('meanNear')
% 
% subplot(plotNum,3,15);
% bar(meanSpin)
% 
%  title('meanSpin')
% 
% subplot(plotNum,3,18);
% bar(multSpin)
% 
% title('multSpin')
% 
% subplot(plotNum,3,21);
% bar(multNear)
% 
% title('multNear')
% 
% 
% 
% %%
% 
% clf
% subplot(4,1,1)
% bar(sumSyn)
% title('synapse number')
% 
% subplot(4,1,2)
% bar(filFrac)
% title('filFrac')
% 
% 
% subplot(4,1,3)
% bar(multSpin)
% title('multSpin')
% 
% 
% subplot(4,1,4)
% bar(multNear)
% title('multNear')
% 
% 


%% Congraph



preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

minCon = 5;
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
    meanInPre(i) = sum(inNum(isPre));
end

for i = 1:length(postNames)
    isPost = postID == postNames(i);
    meanInPost(i)  = sum(inNum(isPost));
end

[sortedPreMean idx] = sort(meanInPre,'descend');
preNames = preNames(idx);
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

for i = 1:length(inCheck)
    subplot(length(inCheck),1,i);
    image(conI(:,:,i)*20)
end


conRat = conI(:,:,2)./conI(:,:,1);
clf
image(conRat*100)

colCon(:,:,1) = conI(:,:,1)*30;
colCon(:,:,2) = sum(conI(:,:,3:end),3)*100;
colCon(:,:,3) = sum(conI(:,:,5:end),3)*0;

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
colCon(:,:,1) = conI(:,:,1)*30;
colCon(:,:,2) = sum(conI(:,:,3:end),3)*100;
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


intersect(isBig, conTo(3).rgcList)
