    MPN = GetMyDir
    load([MPN 'obI.mat'])

    targCells = [108 201];
    conTo = makeConTo(obI,targCells);
    
    
    synProp = obI.nameProps.synProp
    preID = [synProp.pre];
    postID = [synProp.post];
    
    preNames = [1028 1032 2033 2034 2035 2003 2007]
    %postNames = unique(postID);
    
    tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);

    
    %% Filter post
    usePost = postID * 0;
    for i = 1:length(postID)
        
      %usePost(i) =   (sum(tcrList == postID(i)) >0) & (postID(i)>200);
       usePost(i) =   synProp(i).tcr > 0;

    end
    
    postNames = intersect(unique(postID(usePost>0)),[conTo.tcrList]);
    
    
    %% Plot Pres
    plotNum = length(preNames)
    boutBin = 0:10;
    for a = 1:plotNum;
        
        for p = 1:length(postNames)
            
        
        getProps = (preID == preNames(a)) & (postID == postNames(p));
        synNum(a,p) = sum(getProps);
        nearBout(a,p) = mean([synProp(getProps).boutonNeighbor]); 
        spinNum(a,p) = mean([synProp(getProps).spineNum]);
        filFrac(a,p) = mean([synProp(getProps).filimentous]);
        axAxFrac(a,p) = mean([synProp(getProps).axoAxonic]);
        
        end
        
        
        isSyn = synNum(a,:)>0;
        
        
        
        
        %%X
        subplot(plotNum,3,(a-1)*3+1);
        scatter(synNum(a,isSyn),nearBout(a,isSyn));
        %xlim([min(boutBin)-1 max(boutBin)])
        
         subplot(plotNum,3,(a-1)*3+2);
                scatter(synNum(a,isSyn),spinNum(a,isSyn));

         %       xlim([min(boutBin)-1 max(boutBin)])

        
         subplot(plotNum,3,(a-1)*3+3);
         
        scatter(synNum(a,isSyn),filFrac(a,isSyn))
        
    end
    
    %%
    clf
    subplot(3,1,1)
    scatter(synNum(synNum>0),nearBout(synNum>0))
    subplot(3,1,2)
    scatter(synNum(synNum>0),spinNum(synNum>0))
     subplot(3,1,3)
    scatter(synNum(synNum>0),filFrac(synNum>0))
    
    
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
        title('filFrac')
        bar(filFrac)
        
        
         subplot(plotNum,3,6);
         title('inhibFrac')
        bar(inhibFrac)
        
        subplot(plotNum,3,9);
        title('sumSyn')
        bar(sumSyn)
    
        subplot(plotNum,3,12);
        title('meanNear')
        bar(meanNear)
        
        subplot(plotNum,3,15);
        title('meanSpin')
        bar(meanSpin)
        
        subplot(plotNum,3,18);
        title('multSpin')
        bar(multSpin)
        
        subplot(plotNum,3,21);
        title('multNear')
        bar(multNear)
    
    
    
    
    
    
    
    
    
    
    
    

    
    