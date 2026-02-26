


    clear all
    %MPN = GetMyDir;
    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    cellList = unique([plusOne; useList.preList(:); useList.postList(:)]);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
    
    
    %%
    edges = obI.nameProps.edges;    
    anchors = obI.colStruc.anchors(edges(:,3),:);
    
    useEdge = sum(anchors>0,2)>0;
    
    edges = edges(useEdge,:);
    anchors = anchors(useEdge,:);
    
    
    pixSize = obI.em.res .* [4 4 1];
    
    for i = 1:3
        anchors(:,i) = anchors(:,i) .* pixSize(i)/1000;
    end
    
    min(anchors,[],1)
    max(anchors,[],1)
    
    %% Grab synapses that are not on seeds and group them according to seed
    
    badSeed = [];
    for i = 1:length(seedList);
        badSeed = cat(1,badSeed,find(edges(:,1) == seedList(i)));
        
        ax2seed = seedPref.ax2seed(i,:);
        axList = seedPref.axList(ax2seed>0);
        
        axSyn = [];
        for a = 1:length(axList)
            axSyn = cat(1,axSyn,find(edges(:,2)==axList(a)));
        end
        
        
        synGroup{i} = unique(axSyn);
        
    end
    
    for i = 1:length(synGroup)
        synGroup{i} = setdiff(synGroup{i}, badSeed);
    end
    
    
    %% show all
    
    
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
    [samePreHist binpos] = hist(dists(samePre>0),hRange);
    samePostHist = hist(dists(samePost>0),hRange);
    difPreHist = hist(dists(samePre==0),hRange);
    difPostHist = hist(dists(samePost==0),hRange);

    
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
    [samePreHist binpos] = hist(dists(samePre>0),hRange);
    samePostHist = hist(dists((samePost>0) & (samePre==0)),hRange);
    difPreHist = hist(dists(samePre==0),hRange);
    difPostHist = hist(dists((samePost==0) & (samePre==0)),hRange);

    
    clf
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    
   
        bar(hRange,samePostHist./(difPostHist+samePostHist),'b')
    
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
    [samePreHist binpos] = hist(dists(samePre>0),hRange);
    samePostHist = hist(dists(samePost>0),hRange);
    difPreHist = hist(dists(samePre==0),hRange);
    difPostHist = hist(dists(samePost==0),hRange);
    

    
    subplot(length(synGroup),length(synGroup),(g1-1)*length(synGroup) + g2)
    plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
    hold on
    plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
    hold off
    ylim([0 1])
    title(sprintf('groups %d vs %d',seedList(g1), seedList(g2)));
        
        
        
        end
    end
    
    
    
    
    
    
    
    
    
    
    