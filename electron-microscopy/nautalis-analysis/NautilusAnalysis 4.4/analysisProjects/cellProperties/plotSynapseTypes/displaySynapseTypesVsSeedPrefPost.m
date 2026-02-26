
load([MPN 'obI.mat'])

seedList = [108 201];
preNames = [1028 1032 2033 2034 2035 2003 2007]


useList = obI2cellList_seedInput(obI,seedList);
conTo = makeConTo(obI,seedList);
seedPref = seedPreferences(seedList,useList);

synProp = obI.nameProps.synProp
preID = [synProp.pre];
postID = [synProp.post];
uPost = unique(intersect(postID(postID>0),useList.postList))

countPost = histc(postID,uPost);

%% Filter post

tcrPref = seedPref.cellPrefAxExcluded(:,:,1)./sum(seedPref.cellPrefAxExcluded,3);
prefOK = ~isnan(tcrPref);

%% Plot Pref
plotNum = 3;
boutBin = 0:10;

colors = [ 1 0 0; .7 .3 0; 0 1 0; .3 .7 0; 0 .7 .3; 0 0 1; 0 .3 .7]
colormap(colors)
for p = 1:length(uPost);
     filFrac = zeros(length(preNames),1);
     axAxFrac = zeros(length(preNames),1);
     dendFrac = zeros(length(preNames),1);

    for a = 1:length(preNames)
        
        %         axTarg = useList.preList==preNames(a);
        %         postTarg = uPost(p);
        %         useTCR{1} = useList.postList((tcrPref(axTarg,:)==0) & (prefOK(axTarg,:)));
        %         useTCR{2} = useList.postList((tcrPref(axTarg,:)>0) &(tcrPref(axTarg,:)<1) & (prefOK(axTarg,:)));
        %         useTCR{3} = useList.postList((tcrPref(axTarg,:)==1) & (prefOK(axTarg,:)));
        %        
        %         for c = 1:length(colors)
        %
        %             usePost = postID * 0;
        %             for i = 1:length(postID)
        %                 usePost(i) = sum(useTCR{c} == postID(i)) >0;
        %             end
        
        getProps = (preID == preNames(a)) & (postID == uPost(p));
        nearBout(:,a) = hist([synProp(getProps).boutonNeighbor],boutBin);
        spinNum(:,a) = hist([synProp(getProps).spineNum],boutBin);
        filFrac(a) = mean(double([synProp(getProps).filimentous]));
        axAxFrac(a) = mean(double([synProp(getProps).axoAxonic]));
        dendFrac(a) = mean(double([synProp(getProps).isDend]));
        synCount(a) = sum(getProps);

    end
    filFrac(isnan(filFrac)) = 0;
    axAxFrac(isnan(axAxFrac)) = 0;
    dendFrac(isnan(dendFrac)) = 0;
    
    subplot(6,1,1);
    bar(boutBin,nearBout)
    xlim([min(boutBin)-1 max(boutBin)])
    
    subplot(6,1,2);
    bar(boutBin,spinNum)
    xlim([min(boutBin)-1 max(boutBin)])
        
    subplot(6,1,3);
    bar(synCount)

    
    subplot(6,1,4);
    bar(filFrac)
    ylim([0 1])

     
    subplot(6,1,5);
    bar([dendFrac])
    ylim([0 1])

    
    subplot(6,1,6)
    bar([axAxFrac])
    ylim([0 1])
    pause
    
end