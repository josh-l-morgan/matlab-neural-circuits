clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])
load([MPN 'butSize2.mat'])

seedList = [108 201];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
postList = useList.postList;
conTo = makeConTo(obI,seedList);


synPropRaw = obI.nameProps.synProp
clear synProp
c = 0;
for i = 1:length(synPropRaw)
    if ~isempty(synPropRaw(i).pre)
        c = c+1;
        synProp(c) =synPropRaw(i)
    end
end
preID = [synProp.pre];
postID = [synProp.post];
obID = [synProp.objectID]
obI.nameProps.names{obID}

preNames = [1028 1032 2033 2034 2035 2003 2007]
%postNames = [130]; %unique(postID);
%postNames = unique(postID);

tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);


%% Filter post
usePost = postID * 0;
useSeed = usePost;
for i = 1:length(postID)
    
    %usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>200);
    %usePost(i) =   synProp(i).tcr;
%     usePost(i) =  sum(postNames == postID(i));
%     usePost(i) =  sum(postNames == postID(i));
%     usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>1000);
%     usePost(i) = sum(setdiff(postList,seedList) == postID(i));
%         useSeed(i) = sum(seedList == postID(i));
    usePost(i) = sum(postList == postID(i));

end
%usePost = postID * 0+1;


%% Plot Pres
plotNum = length(preNames)
boutBin = 0:10;
axList = butSize.axList;

for a = 1:plotNum;
    
    getProps = (preID == preNames(a)) & (usePost>0);
    nearBout = hist([synProp(getProps).boutonNeighbor],boutBin);
    spinNum = hist([synProp(getProps).spineNum],boutBin);
    filFrac = mean([synProp(getProps).filimentous]);
    axAxFrac = mean([synProp(getProps).axoAxonic]);
    
    subplot(plotNum,3,(a-1)*3+1);
    bar(boutBin,nearBout)
    xlim([min(boutBin)-1 max(boutBin)])
    
    subplot(plotNum,3,(a-1)*3+2);
    bar(boutBin,spinNum)
    xlim([min(boutBin)-1 max(boutBin)])
    
    
    subplot(plotNum,3,(a-1)*3+3);
    
    bar([filFrac 0])
    hold on
    bar([0 axAxFrac],'r')
    hold off
    ylim([0 1])
    
end



%% Plot Pres
clf
plotNum = length(preNames)
boutBin = 0:10;
pcol = {'r','g','b','c','m','y','k'}
for a = 1:plotNum;
    
    axTarg = find(axList == preNames(a));
    butVols = butSize.butVols{axTarg};
    butDiam = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2;
    meanDiameter(a) = mean(butDiam);
    SE(a) = std(butDiam)/sqrt(length(butDiam));
    
    getProps = (preID == preNames(a)) & usePost;
    plot(getProps,pcol{a})
    hold on
    
    filFrac(a) = mean([synProp(getProps).filimentous]);
    sumSyn(a) = sum(getProps);
    %axAxFrac = mean([synProp(getProps).axoAxonic]);
    
    meanNear(a) = mean([synProp(getProps).boutonNeighbor]);
    stdNear(a) = std([synProp(getProps).boutonNeighbor]);%/sqrt(sum(getProps));
    meanSpin(a) = mean([synProp(getProps).spineNum]);
    stdSpin(a) = std([synProp(getProps).spineNum]);%/sqrt(sum(getProps));
    
    multSpin(a) = mean([synProp(getProps).spineNum]>1);
    multNear(a) =  mean([synProp(getProps).boutonNeighbor]>1);
    
end
hold off

[meanDiameter' SE']

clf
subplot(5,1,1)
bar(sumSyn,'k')
title('synapse number')

subplot(5,1,2)
bar(meanDiameter,'k')
title('filFrac')

subplot(5,1,3)
bar(filFrac,'k')
title('filFrac')


subplot(5,1,4)
bar(multSpin,'k')
title('multSpin')


subplot(5,1,5)
bar(multNear,'k')
title('multNear')












