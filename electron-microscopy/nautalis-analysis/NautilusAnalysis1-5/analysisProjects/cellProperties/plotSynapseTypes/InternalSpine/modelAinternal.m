
     
     %% Monte carlo test of frequency of internal spine formation on cell A
     %%Data comes from ...
     
          inList = ax.inList;
          postList = ax.postIDs;
          
          targCell = 108;
      preHist = preTo(allEdges,targCell);
     useAx = intersect(ax.axID,preHist(:,1)
     
     for i = 1:length(useAx)
         
         targ = find(ax.axID == useAx(i));
         allSyn(i) = ax.synNum(targ);
         anyNum(i)  =ax.anyNum(targ);
         inAx = ax.inList{targ};
         postAx = ax.postIDs{targ};
         
         offInt{i} = inAx(~(postAx==targCell));
         useOff(i) = sum(~(postAx==targCell))>=3;
         onInt{i} = inAx((postAx==targCell));

         targ = find(preHist(:,1)==useAx(i));
         onSyn(i) = preHist(targ,2);
         
     end
     
     useAx = find(useOff);
     
     realInt = 0;
     for i = 1:length(useAx)
        realInt = realInt + sum(onInt{useAx(i)})
     end
     
     
    reps = 10000; 
    isPerf = zeros(reps,length(useAx));
     for r = 1:reps
        for a = 1:length(useAx)
            targ = useAx(a);
            pickRand =  ceil(rand(onSyn(targ),1) * length(offInt{targ}));
            gotSyn = offInt{targ}(pickRand);
           isPerf(r,a) = sum(gotSyn); 
        end
     end
     totPerf = sum(isPerf,2);
     
     P = sum(totPerf<=realInt)/length(totPerf)
     median(totPerf)
     rangeX(totPerf)
         
     
     
     %%
     %{
     
minSyn = 1;
ratHistBin = [0:.1:1];
clear percentGiant
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

%%
    preHist = preTo(allEdges,108);
     
     useAx = intersect(ax.axID,preHist(:,1))
     for i = 1:length(useAx)
         
         targ = find(ax.axID == useAx(i));
         allSyn(i) = ax.synNum(targ);
         anyNum(i)  =ax.anyNum(targ);

         targ = find(preHist(:,1)==useAx(i));
         onSyn(i) = preHist(targ,2);
         
     end
     offSyn = allSyn-onSyn;
     
     %%
     
     reps = 10000;
     isPerf = zeros(reps,length(onSyn));
     for r = 1:reps
        for a = 1:length(onSyn)
           pickRand =  randperm(allSyn(a),onSyn(a));
           isPerf(r,a) = sum(pickRand<=anyNum(a)); 
        end
     end
     totPerf = sum(isPerf,2);
     
     P = sum(totPerf<=1)/length(totPerf)
     median(totPerf)
     rangeX(totPerf)
     
     %%
     
     useSyn = (offSyn>0)
     useOff = offSyn(useSyn);
     useAny = anyNum(useSyn);
     useOn = onSyn(useSyn);
     anyFrac = useAny./useOff;
     
     
     reps = 10000;
     isPerf = zeros(reps,length(onSyn));
     for r = 1:reps
        for a = 1:length(useOff)
           isPerf(r,a) = sum(rand(useOn(a),1)<=anyFrac(a)); 
        end
     end
     totPerf = sum(isPerf,2);
     
     P = sum(totPerf<=1)/length(totPerf)
     median(totPerf)
     rangeX(totPerf)
     
     
     %}
     
         
         
     
     
     