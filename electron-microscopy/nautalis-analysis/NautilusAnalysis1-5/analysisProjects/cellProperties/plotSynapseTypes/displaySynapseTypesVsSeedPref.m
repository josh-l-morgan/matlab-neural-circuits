    clear all
    
    MPN = GetMyDir
    load([MPN 'obI.mat'])

    seedList = [108 201 170];
    preNames = [1028 1032 2033 2034 2035 2003 2007]

   
    useList = obI2cellList_seedInput(obI,seedList);
    conTo = makeConTo(obI,seedList);
    seedPref = seedPreferences(seedList,useList);
    
    synProp = obI.nameProps.synProp
    preID = [synProp.pre];
    postID = [synProp.post];
    
    %% Filter post
   
    
    tcrPref = seedPref.cellPrefAxExcluded(:,:,1)./sum(seedPref.cellPrefAxExcluded,3);
    prefOK = ~isnan(tcrPref);
    
    
    
    %% Plot Pref
    plotNum = length(preNames)
    boutBin = 0:10;
    
    colors = {'r','g','b'}
    prefThresh = [10 10];
    for a = 1:plotNum;
        
        axTarg = useList.preList==preNames(a);
        
        useTCR{1} = useList.postList((tcrPref(axTarg,:)<=prefThresh(1)) & (prefOK(axTarg,:)));
        useTCR{2} = useList.postList((tcrPref(axTarg,:)>prefThresh(1)) &(tcrPref(axTarg,:)<prefThresh(2)) & (prefOK(axTarg,:)));
        useTCR{3} = useList.postList((tcrPref(axTarg,:)>=prefThresh(2)) & (prefOK(axTarg,:)));
        filFrac = zeros(10,3);
        axAxFrac = zeros(10,3);
        for c = 1:length(colors)
           
            usePost = postID * 0;
            for i = 1:length(postID)
                usePost(i) = sum(useTCR{c} == postID(i)) >0;
            end
            
        getProps = (preID == preNames(a)) & usePost;
        nearBout(:,c) = hist([synProp(getProps).boutonNeighbor],boutBin); 
        spinNum(:,c) = hist([synProp(getProps).spineNum],boutBin);
        filFrac(1,c) = mean([synProp(getProps).filimentous]);
        axAxFrac(2,c) = mean([synProp(getProps).axoAxonic]);
        synNum = sum(getProps);
        end
        filFrac(isnan(filFrac)) = 0;
        axAxFrac(isnan(axAxFrac)) = 0;
        
        subplot(plotNum,3,(a-1)*3+1);
        bar(boutBin,nearBout)
        xlim([min(boutBin)-1 max(boutBin)])
        
         subplot(plotNum,3,(a-1)*3+2);
        bar(boutBin,spinNum)
                xlim([min(boutBin)-1 max(boutBin)])

        
         subplot(plotNum,3,(a-1)*3+3);
         
        bar(filFrac)
        hold on
        bar([axAxFrac],'r')
        hold off
        ylim([0 1])
        xlim([.5 2.5])
        
        
    end