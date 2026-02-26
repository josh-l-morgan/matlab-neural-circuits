clear all
load('MPN.mat')
load([MPN 'butSize2.mat'])
load([MPN 'obI.mat'])

%%
seedList = [108 201 109 903 907];
seedNum = length(seedList);
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

postList = useList.postList;
edges = butSize.edges;
con = edge2con(edges);
axList = butSize.axList;

conPref =seedPreferences(seedList,useList);
isMixed = sum(conPref.sharedSyn([1 2],:)>0,1)==2;
mixList = conPref.cellList(isMixed);
%mixList = postList(randperm(length(postList),26));
unMixList = setdiff(postList,mixList);

%%
clear compSizes
for a = 1:length(axList)
   compSizes(a,:) = cat(2,a, size(butSize.butVols{a},2),...
   size(butSize.axProfile{a},1),...
   size(butSize.synDat(a).edges,1));
end

clear preSeed onSeedVols offSeedVols onMixVols
for s = 1:seedNum % check each seed
    preSeed{s} = unique(edges(edges(:,1) == seedList(s),2));
    preSeed{s} = intersect(preSeed{s},axList);
    onSeedVols{s} = [];
    offSeedVols{s} = [];
    onMixVols{s} = [];
    for a = 1:length(preSeed{s}) % check each axon
        
        axID = preSeed{s}(a);
        axTarg = find(axList == axID);
        butVols = butSize.butVols{axTarg};
        
        butVols = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2
        %butVols = butVols * butSize.voxVol;
        
        postBut = butSize.synDat(axTarg).edges(:,1);
        onSeed = find(postBut == seedList(s));
        offSeed = [];
        for p = 1:length(unMixList)
            if (unMixList(p) ~=seedList(s))
                offSeed = [offSeed; find(postBut == unMixList(p))];
            end
        end
        onMix = [];
        for p = 1:length(mixList)
            if (mixList(p) ~=seedList(s))
                onMix = [onMix; find(postBut == mixList(p))];
            end
        end
        onMixVols{s} = [onMixVols{s} butVols(onMix)];
        onSeedVols{s} = [onSeedVols{s} butVols(onSeed)];
        offSeedVols{s} = [offSeedVols{s} butVols(offSeed)];
    end
    onMixVols{s} = onMixVols{s}(onMixVols{s}>0);
    onSeedVols{s} = onSeedVols{s}(onSeedVols{s}>0);
    offSeedVols{s} = offSeedVols{s}(offSeedVols{s}>0);
    
end

%%
allBut = [offSeedVols{:} onSeedVols{:} onMixVols{:}];%[butSize.butVols{:}] * butSize.voxVol;
%butVolHistBin = [1:2:max(allBut)];
%butVolHistBin = [0:2:max(allBut)];
binWidth = .3;
butVolHistBin = [binWidth/2:binWidth:max(allBut)+binWidth];

xmax = max(allBut) + binWidth ;
plotSeed = [1 2];
clear collectVols 
collectBars = []

for s = 1:length(plotSeed)
%     
%     totBut = length(onMixVols{s}) + length(onSeedVols{s}) + length(offSeedVols{s})
%     length(postBut)
    
%%Plot on
    subplot(length(plotSeed),3,s*3-2)
    foundSyn = find(onSeedVols{s});
    onSyn = onSeedVols{s}(foundSyn);
    histVol = hist(onSyn,butVolHistBin);
    bar(butVolHistBin,histVol,'k')    
    %ylim([0 60])
    xlim([-binWidth xmax])
    collectVols{s,1} = onSyn;
    collectBars{1}(s,:) = histVol;
    N = sum(histVol)
    
    
    
    
    
    %%plot off Seed
    subplot(length(plotSeed),3,s*3-1);
    foundSyn = find(offSeedVols{s});
    offSyn = offSeedVols{s}(foundSyn);
    histVol = hist(offSyn,butVolHistBin);
    bar(butVolHistBin,histVol,'k')  
    %ylim([0 60])
          xlim([-binWidth xmax])
    collectBars{2}(s,:) = histVol;
    N = sum(histVol)

      
    %%plot mixed
    subplot(length(plotSeed),3,s*3);
    foundSyn = find(onMixVols{s});
    mixSyn = offSeedVols{s}(foundSyn);
    histVol = hist(mixSyn,butVolHistBin);
    bar(butVolHistBin,histVol,'k')  
    %ylim([0 60])
          xlim([-binWidth xmax])
    collectBars{3}(s,:) = histVol;
    N = sum(histVol)

    %%Stats on
    meanOn = mean(onSyn)
    Non = length(onSyn)
    range95on = rangeX(onSyn,0.95)
    SEon = std(onSyn)/sqrt(length(onSyn))
    
    %%Stats off
    meanOff = mean(offSyn)
    Noff = length(offSyn)
    range95off = rangeX(offSyn,0.95)
    SEoff = std(offSyn)/sqrt(length(offSyn))

    %%Stats mixed
    meanMixed = mean(mixSyn)
    Nmix = length(mixSyn)
    range95mix = rangeX(mixSyn,0.95)
    SEmix = std(mixSyn)/sqrt(length(mixSyn))
    
    
    collectVols{s,2} = offSyn;
    collectVols{s,3} = mixSyn;

end
    

ranksum(collectVols{1,1},collectVols{2,1})
ranksum(collectVols{1,2},collectVols{2,2})
ranksum(collectVols{1,3},collectVols{2,3})

pause(.5)
%% merge bars

for b = 1:3
    subplot(3,1,b)
    bar(butVolHistBin,collectBars{b}','barwidth',1.6)
end

%% plot linked tcr from seed perspective


clf

for s = 1:length(seedList)

    preSeed = intersect(axList,conTo(s).rgcList);
    linkedTC = intersect(postList,conTo(s).tcrList);
    
    postSizes = {};
    for a = 1:length(preSeed)
        axTarg = find(butSize.axList == preSeed(a))
        butVol = butSize.butVols{axTarg};
        butVol = (butVol*3/4/pi).^(1/3)*butSize.voxLength * 2;
        
        postBut = butSize.synDat(axTarg).edges(:,1);
        for p = 1:length(postBut)
           targ =  find(linkedTC == postBut(p));
           if ~isempty(targ)
            try
               postSizes{targ} = cat(2,postSizes{targ},butVol(p));
            catch err
                postSizes{targ} = butVol(p);
            end
           end
        end
    end
    clear meanTC totTC
    for t = 1:length(postSizes)
       meanTC(t) = mean(postSizes{t}); 
       totTC(t) = length(postSizes{t});
    end
    
    
    subplot(length(seedList),1,s)
    scatter(totTC,meanTC,20,'k','filled');
    hold on
    targ = find(linkedTC == seedList(s));
    scatter(totTC(targ),meanTC(targ),40,'r','filled');
    plot([0 230],[meanTC(targ) meanTC(targ)],'r')
    %ylim([0 4000])
    xlim([0 230])
    

end

% 




%% Make get list table of diameters
clear meanDiam
for a = 1:length(butSize.axList);
    butVols = butSize.butVols{a};
    butDiam = vol2diam(butVols * butSize.voxVol);
    meanDiam(a) = mean(butDiam);
end
propMat = [butSize.axList' meanDiam']

%% Make get list table of diameters seed not excluded
clear meanDiam
sumPost = zeros(1,length(postList));
countPost = zeros(1,length(postList));
for a = 1:length(butSize.axList);
    synDat = butSize.synDat(a);
    postSyn = synDat.edges(:,1);
    butVols = butSize.butVols{a};
    butDiam = vol2diam(butVols * butSize.voxVol);
    sumAx = 0;
    countAx = 0;
    for p = 1:length(postSyn)
       targ = find(postList == postSyn(p));
       if ~isempty(targ)
          sumPost(targ) = sumPost(targ) + butDiam(p);
          countPost(targ) = countPost(targ) + 1;
       end
       
       if 1%sum(seedList == postSyn(p)) == 0
           sumAx = sumAx + butDiam(p);
           countAx = countAx + 1;
       end
        
    end
   
    meanDiam(a) = sumAx/countAx;
    axSynCount(a) = countAx;
end
propMat = [butSize.axList' meanDiam' axSynCount';
    postList' sumPost'./countPost' countPost']



%% show axons
clf
axNum = length(butSize.axList);
 for a = 1:axNum % check each axon
         butVols = butSize.butVols{a};
        
        butVols = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2
        %butVols = butVols * butSize.voxVol;
        postBut = butSize.synDat(axTarg).edges(:,1);
        
        
        subplot(3,3,mod(a,9)+1)
    histVol = hist(butVols,butVolHistBin);
    bar(butVolHistBin,histVol,'k')    
    %ylim([0 60])
    xlim([-binWidth 6])
    N = sum(histVol)
    pause(.01)
        
        
        
        
    end





