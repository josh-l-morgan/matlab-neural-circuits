


    clear all
    %MPN = GetMyDir;
    %load('MPN.mat')
    MPN = 'C:\Users\jlmorgan\Documents\Segmentations\AndermannChen\export\MAC_merge_mat\';

    load([MPN 'obI.mat']);
   
%     useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
%     cellList = unique([plusOne; useList.preList(:); useList.postList(:)]);
%     
%     seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
    
    
    %% get synapse positions
    edges = obI.nameProps.edges;    
    anchors = obI.colStruc.anchors(edges(:,3),:);
    
    useEdge = sum(anchors>0,2)>0;
    useEdge = useEdge & ~sum(edges == 0,2);
    
    edges = edges(useEdge,:);
    anchors = anchors(useEdge,:);
    numSyn = size(edges,1);

    emRes = [6 4 30]; %obI.em.res;
    pixSize = emRes .* [4 4 1];
    
    for i = 1:3
        anchors(:,i) = anchors(:,i) .* pixSize(i)/1000;
    end
    
    min(anchors,[],1)
    max(anchors,[],1)
    
       
      clf

    
 %% synapse distance matrix
 
    %%Get all synapses (not on seeds)   
    syn1 = 1:numSyn;
    syn2 =  1:numSyn;
    
    
    subplot(3,1,1)
    
    
    colmap = hsv(100);
    postCells = unique(edges(:,1));
    %postCells = postCells(postCells>0);
    numPost = length(postCells);
    pickPost = round((1:numPost) * 100/numPost);
    pickPost = pickPost(randperm(numPost));
    cellCol = colmap(pickPost,:);
    lookupCol(postCells,:) = cellCol;
    
    memCol.postCells = postCells;
    memCol.cellCol = cellCol;
    memCol.lookupCol = lookupCol;
    
%     set(gcf,colormap(lookupCol))
%     scatter(anchors(syn1,1),anchors(syn1,2),50,edges(syn1,1),'o','filled')
%subplot(1,1,1)
    for i = 1:length(syn1)
        scatter(anchors(syn1(i),1),anchors(syn1(i),2),45,lookupCol(edges(syn1(i),1),:),...
            'filled','markeredgecolor','w')
        %set(gca,'MarkerEdgeColor',lookupCol(edges(syn1(i))));
        %scatter(anchors(syn1(i),1),anchors(syn1(i),2),lookupCol(edges(syn1(i),1),:))
        %scatter(anchors(syn1(i),1),anchors(syn1(i),2),50,edges(1,1),'o','filled')
        hold on
    end
    hold off
    xlim([140 210])
    ylim([220 290])
    
    samePost = zeros(length(syn1),length(syn2));
    samePre = samePost; dists = samePost;useComp = samePost;
    for s1 = 1:length(syn1)
        for s2 = 1:length(syn2)
            
            samePost(s1,s2) = edges(syn1(s1),1)==edges(syn2(s2),1);
            samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
            dists(s1,s2) = sqrt((anchors(syn1(s1),1) - anchors(syn2(s2),1)).^2 + ...
                (anchors(syn1(s1),2) - anchors(syn2(s2),2)).^2 + ...
                (anchors(syn1(s1),3) - anchors(syn2(s2),3)).^2);
            useComp(s1,s2) = s1<s2;
            sameComp(s1,s2) = s1==s2;
       end
    end
    
    symDists = dists;
    symDists(sameComp) = inf;
    dists(~useComp) = inf; 
    
    hRange = [0:6:100];
    showRange = (hRange(2:end)+hRange(1:end-1))/2;
        %showRange = hRange(1:end-1);

    [samePreHist] = histcounts(dists(samePre>0),hRange);
    samePostHist = histcounts(dists(samePost>0),hRange);
    difPreHist = histcounts(dists(samePre==0),hRange);
    difPostHist = histcounts(dists(samePost==0),hRange);

    
    subplot(3,1,2)
    sumBoth = sum(samePostHist) + sum(difPostHist);
    bar(showRange(1:end-1),cat(1,difPostHist(1:end-1)/sumBoth,...
        samePostHist(1:end-1)/sumBoth)', 1.5)
    xlim([-2 100])
    hold on
    subplot(3,1,3)
    %plot(hRange(1:end-1),samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1)),'b')
    bar(showRange(1:end-1),samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1)),'b')

    xlim([-2 100])
    ylim([0 1])
    hold off
    ylim([0 1])
    
   
       % bar(hRange,samePostHist./(difPostHist+samePostHist),'b')
    
    
return        
 %% Find clustered synapses
 
 synWithin = 10; %threshold for clustering
 symDists = dists + dists';
 
 
 synThresh = (dists<=synWithin);
 sameThresh = synThresh .* samePre;
 
 image(sameThresh * 1000)
 
 isClustered = sum(sameThresh,1) + sum(sameThresh,2)';
 histcounts(isClustered,[0:1:10])
 
 
 
   hRange = [0:1:100];
    [samePreHist binpos] = histcounts(dists((samePre>0)&(sameThresh)),hRange);
    samePostHist = histcounts(dists((samePost>0)&(sameThresh)),hRange);
    difPreHist = histcounts(dists((samePre==0)&(sameThresh)),hRange);
    difPostHist = histcounts(dists((samePost==0)&(sameThresh)),hRange);

    
    clf
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    
   
        bar(hRange,samePostHist./(difPostHist+samePostHist),'b')
 
 
    

      %% show all, diff axon
    
    
    syn1 = unique(cat(1,synGroup{:}));
    syn2 =  unique(cat(1,synGroup{:}));
    
    
    samePost = zeros(length(syn1),length(syn2));
    samePre = samePost; dists = samePost;useComp = samePost;
    for s1 = 1:length(syn1)
        for s2 = 1:length(syn2)
            
            samePost(s1,s2) = edges(syn1(s1),1)==edges(syn2(s2),1);
            samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
            dists(s1,s2) = sqrt((anchors(syn1(s1),1) - anchors(syn2(s2),1)).^2 + ...
                (anchors(syn1(s1),2) - anchors(syn2(s2),2)).^2 + ...
                (anchors(syn1(s1),3) - anchors(syn2(s2),3)).^2);
            useComp(s1,s2) = s1<s2;
       end
    end
    
    dists(~useComp) = inf;
    
    hRange = [0:1:100];
    [samePreHist binpos] = histcounts(dists(samePre>0),hRange);
    samePostHist = histcounts(dists((samePost>0) & (samePre==0)),hRange);
    difPreHist = histcounts(dists(samePre==0),hRange);
    difPostHist = histcounts(dists((samePost==0) & (samePre==0)),hRange);

    
    clf
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    
   
       % bar(hRange,samePostHist./(difPostHist+samePostHist),'b')
    
    %% compare groups
    
    for g1 = 1:length(synGroup);
        for g2 = 1:length(synGroup);
        
        
        
    
    syn1 = unique(cat(1,synGroup{g1}));
    syn2 =  unique(cat(1,synGroup{g2}));
    
    
    samePost = zeros(length(syn1),length(syn2));
    samePre = samePost; dists = samePost; useComp = samePost;
    for s1 = 1:length(syn1)
        for s2 = 1:length(syn2)
            
            samePost(s1,s2) = edges(syn1(s1),1)==edges(syn2(s2),1);
            samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
            dists(s1,s2) = sqrt((anchors(syn1(s1),1) - anchors(syn2(s2),1)).^2 + ...
                (anchors(syn1(s1),2) - anchors(syn2(s2),2)).^2 + ...
                (anchors(syn1(s1),3) - anchors(syn2(s2),3)).^2);
            useComp(s1,s2) = s1<s2;
       end
    end
    
    %dists(~useComp) = inf;
    sameBoth = samePre & samePost;
    
    hRange = [0:1:100];
    [samePreHist binpos] = histcounts(dists(samePre>0),hRange);
    samePostHist = histcounts(dists(samePost>0),hRange);
    difPreHist = histcounts(dists(samePre==0),hRange);
    difPostHist = histcounts(dists(samePost==0),hRange);
    

    
    subplot(length(synGroup),length(synGroup),(g1-1)*length(synGroup) + g2)
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    title(sprintf('groups %d vs %d',seedList(g1), seedList(g2)));
        
        
        
        end
    end
    
    
     %% compare AA, AB, BB
    
    for g1 = 1:length(synGroup);
        for g2 = 1:length(synGroup);
        
        
        
    
    syn1 = unique(cat(1,synGroup{g1}));
    syn2 =  unique(cat(1,synGroup{g2}));
    
    
    samePost = zeros(length(syn1),length(syn2));
    samePre = samePost; dists = samePost; useComp = samePost;
    for s1 = 1:length(syn1)
        for s2 = 1:length(syn2)
            
            samePost(s1,s2) = edges(syn1(s1),1)==edges(syn2(s2),1);
            samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
            dists(s1,s2) = sqrt((anchors(syn1(s1),1) - anchors(syn2(s2),1)).^2 + ...
                (anchors(syn1(s1),2) - anchors(syn2(s2),2)).^2 + ...
                (anchors(syn1(s1),3) - anchors(syn2(s2),3)).^2);
            useComp(s1,s2) = s1<s2;
       end
    end
    
    %dists(~useComp) = inf;
    sameBoth = samePre & samePost;
    
    hRange = [0:20:100];
    [samePreHist binpos] = histcounts(dists(samePre>0),hRange);
    samePostHist = histcounts(dists(samePost>0),hRange);
    difPreHist = histcounts(dists(samePre==0),hRange);
    difPostHist = histcounts(dists(samePost==0),hRange);
    

    
    subplot(length(synGroup),length(synGroup),(g1-1)*length(synGroup) + g2)
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    title(sprintf('groups %d vs %d',seedList(g1), seedList(g2)));
        
        
        
        end
    end
    
    
    
    
    
    
    
    
    