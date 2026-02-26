clear all
%MPN = GetMyDir
load('MPN.mat')
load([MPN 'obI.mat'])

targCells = [108 201];
conTo = makeConTo(obI,targCells);
seedList = [108 201 109 903 907]
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
postList = useList.postList;

synProp = obI.nameProps.synProp
preID = [synProp.pre];
postID = [synProp.post];

preNames = [1028 1032 2033 2034 2035 2003 2007]
postNames = [130]; %unique(postID);
%postNames = unique(postID);

tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);


%% Filter post
usePost = postID * 0;
for i = 1:length(postID)
    
    %usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>200);
    %usePost(i) =   synProp(i).tcr;
    usePost(i) =  sum(postList == postID(i));
end


%% Plot Pres
plotNum = length(preNames)
boutBin = 0:10;
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
plotNum = length(preNames)
boutBin = 0:10;
for a = 1:plotNum;
    
    getProps = (preID == preNames(a)) & usePost;
    getAllPre = (preID == preNames(a));
    linNum = sum([synProp(getAllPre).lin]);
    tcrNum = sum([synProp(getAllPre).tcr]);
    inhibFrac(a) = linNum/(linNum+tcrNum);
    
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



subplot(plotNum,3,3);
bar(filFrac)

% title('filFrac')


subplot(plotNum,3,6);
bar(inhibFrac)

% title('inhibFrac')

subplot(plotNum,3,9);
bar(sumSyn)

 title('sumSyn')

subplot(plotNum,3,12);
bar(meanNear)

 title('meanNear')

subplot(plotNum,3,15);
bar(meanSpin)

 title('meanSpin')

subplot(plotNum,3,18);
bar(multSpin)

title('multSpin')

subplot(plotNum,3,21);
bar(multNear)

title('multNear')



%%

clf
subplot(4,1,1)
bar(sumSyn)
title('synapse number')

subplot(4,1,2)
bar(filFrac)
title('filFrac')


subplot(4,1,3)
bar(multSpin)
title('multSpin')


subplot(4,1,4)
bar(multNear)
title('multNear')













