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

conPref = seedPreferences([108 201],useList);
sumSeed = sum(conPref.sharedAx>0,1);
unMix = conPref.cellList(sumSeed == 1);
mix = conPref.cellList(sumSeed >1);



clear usePost useMix useUnMix seedGroup
for i = 1:length(postID)
    
    %usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>200);
    %usePost(i) =   synProp(i).tcr;
%     usePost(i) =  sum(postNames == postID(i));
%     usePost(i) =  sum(postNames == postID(i));
%     usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>1000);
%     usePost(i) = sum(setdiff(postList,seedList) == postID(i));
%         useSeed(i) = sum(seedList == postID(i));
    usePost(i) = sum(postList == postID(i));
    useMix(i) = sum(mix == postID(i));
    useUnMix(i) = sum(unMix == postID(i));
    
    seedGroup(i) = 2-sum([1028 1032] ==preID(i));
    
end
%usePost = postID * 0+1;

mixGroup{2} = useMix;
mixGroup{1} = useUnMix;

axGroup = zeros(length(butSize.synDat),1);
for a = 1:length(butSize.synDat);
    tcList = [butSize.synDat(a).edges(:,1)];
    rgcList = [butSize.synDat(a).edges(:,1)];
    if sum([1028 1032] == butSize.axList(a))
        axGroup(a) = 1;
    elseif sum([ 2033 2034 2035 2003 2007] == butSize.axList(a));
        axGroup(a) = 2;
    end
    clear useMixBut useUnMixBut 

    for i = 1:length(tcList)
        useMixBut(i) = sum(mix == tcList(i));
        useUnMixBut(i) = sum(unMix == tcList(i));
        
    end
    mixGroupBut{a,2} = useMixBut;
    mixGroupBut{a,1} = useUnMixBut;
end

%% Plot Pres
plotNum = length(preNames)
boutBin = 0:10;
axList = butSize.axList;

for a = 1:plotNum;
    
    getProps = (preID == preNames(a)) ;
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
for a = 1:2;
    for m = 1:2
    axTargs = find((axGroup == a)) ;
    useVol = [];
    for t = 1:length(axTargs)
        butVols = butSize.butVols{axTargs(t)};
        butVols = butVols(mixGroupBut{axTargs(t),m}>0);
        useVol = [useVol butVols];
    end
    butDiam = (useVol*3/4/pi).^(1/3)*butSize.voxLength * 2;
    meanDiameter(a,m) = mean(butDiam);
    SE(a,m) = std(butDiam)/sqrt(length(butDiam));
    
    getProps = (seedGroup == a) &  (mixGroup{m}>0);
    plot(getProps,pcol{a})
    hold on
    
    filFrac(a,m) = mean([synProp(getProps).filimentous]);
    sumSyn(a,m) = sum(getProps);
    %axAxFrac = mean([synProp(getProps).axoAxonic]);
    
    meanNear(a,m) = mean([synProp(getProps).boutonNeighbor]);
    stdNear(a,m) = std([synProp(getProps).boutonNeighbor]);%/sqrt(sum(getProps));
    meanSpin(a,m) = mean([synProp(getProps).spineNum]);
    stdSpin(a,m) = std([synProp(getProps).spineNum]);%/sqrt(sum(getProps));
    
    multSpin(a,m) = mean([synProp(getProps).spineNum]>1);
    multNear(a,m) =  mean([synProp(getProps).boutonNeighbor]>1);
    end
end
hold off

[meanDiameter' SE']

clf
subplot(5,1,1)
bar(sumSyn,'k')
title('synapse number')

subplot(5,1,2)
bar(meanDiameter,'k')
title('meanDiameter')

subplot(5,1,3)
bar(filFrac,'k')
title('filFrac')


subplot(5,1,4)
bar(multSpin,'k')
title('multSpin')


subplot(5,1,5)
bar(multNear,'k')
title('multNear')

return

%% Plot Pres
clf
plotNum = length(preNames)
boutBin = 0:10;
pcol = {'r','g','b','c','m','y','k'}
for a = 1:2;
    for m = 1:2
    axTargs = find((axGroup == a)) ;
    useVol = [];
    for t = 1:length(axTargs)
        butVols = butSize.butVols{axTargs(t)};
        butVols = butVols(mixGroupBut{axTargs(t),m}>0);
        useVol = [useVol butVols];
    end
    butDiam = (useVol*3/4/pi).^(1/3)*butSize.voxLength * 2;
    meanDiameter(a,m) = mean(butDiam);
    SE(a,m) = std(butDiam)/sqrt(length(butDiam));
    
    getProps = (seedGroup == a) &  (mixGroup{m}>0);
    plot(getProps,pcol{a})
    hold on
    
    filFrac(a,m) = mean([synProp(getProps).filimentous]);
    poolFil{a,m} = [synProp(getProps).filimentous];
    
    sumSyn(a,m) = sum(getProps);
    %axAxFrac = mean([synProp(getProps).axoAxonic]);
    
    meanNear(a,m) = mean([synProp(getProps).boutonNeighbor]);
    stdNear(a,m) = std([synProp(getProps).boutonNeighbor]);%/sqrt(sum(getProps));
    meanSpin(a,m) = mean([synProp(getProps).spineNum]);
    
    stdSpin(a,m) = std([synProp(getProps).spineNum]);%/sqrt(sum(getProps));
    
    multSpin(a,m) = mean([synProp(getProps).spineNum]>1);
    poolSpin{a,m} = [synProp(getProps).spineNum]>1;

    multNear(a,m) =  mean([synProp(getProps).boutonNeighbor]>1);
    poolNear{a,m} = [synProp(getProps).boutonNeighbor]>1;
    end
end
hold off

[meanDiameter' SE']

clf
subplot(5,1,1)
bar(sumSyn,'k')
title('synapse number')

subplot(5,1,2)
bar(meanDiameter,'k')
title('meanDiameter')

subplot(5,1,3)
bar(filFrac,'k')
title('filFrac')


subplot(5,1,4)
bar(multSpin,'k')
title('multSpin')


subplot(5,1,5)
bar(multNear,'k')
title('multNear')









