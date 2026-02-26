


    clear all
    %MPN = GetMyDir;
    %MPN = 'C:\Users\jlmorgan\Documents\Data\Segmentations\AndermannChen\export\MAC_merge_mat\';

    load('MPN.mat')
    load([MPN 'obI.mat']);
    seedList = [ 108 201 907 903];
    plusOne = 125;
    %seedList = [ 108  201 109 ];
    
    useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
    cellList = unique([plusOne; useList.preList(:); useList.postList(:)]);
    
    seedPref = seedPreferences(seedList,useList);
    allEdges = obI.nameProps.edges(:,[2 1]);
    
    
    
    %% get synapse positions
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
        
        %% found synapses belonging to seed
        badSeed = cat(1,badSeed,find(edges(:,1) == seedList(i)));
        
        %% find axons innervating seed
        ax2seed = seedPref.ax2seed(i,:);
        axList = seedPref.axList(ax2seed>0);
        
        %% find synapses belongin to seeds axons
        axSyn = [];
        for a = 1:length(axList)
            axSyn = cat(1,axSyn,find(edges(:,2)==axList(a)));
        end
        
        
        synGroup{i} = unique(axSyn);
        
    end
    
    %% remove all synapses formed onto seed cells
    for i = 1:length(synGroup)
        synGroup{i} = setdiff(synGroup{i}, badSeed);
    end
    
  
    
 %% synapse distance matrix
 
    %%Get all synapses (not on seeds)   
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
            sameComp(s1,s2) = s1==s2;
       end
    end
    
    symDists = dists;
    symDists(sameComp) = inf;
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
    
    
        
 %% Find clustered synapses
 
 synWithin = 10; %threshold for clustering
 symDists = dists + dists';
 
 
 synThresh = (dists<=synWithin);
 sameThresh = synThresh .* samePre;
 
 image(sameThresh * 1000)
 
 isClustered = sum(sameThresh,1) + sum(sameThresh,2)';
 hist(isClustered,[0:1:10])
 
 
 
   hRange = [0:1:100];
    [samePreHist binpos] = hist(dists((samePre>0)&(sameThresh)),hRange);
    samePostHist = hist(dists((samePost>0)&(sameThresh)),hRange);
    difPreHist = hist(dists((samePre==0)&(sameThresh)),hRange);
    difPostHist = hist(dists((samePost==0)&(sameThresh)),hRange);

    
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
    
    hRange = [0:6:100];
    showRange = (hRange(2:end)+hRange(1:end-1))/2;
    [samePreHist binpos] = histcounts(dists(samePre>0),hRange);
    samePostDifHist = histcounts(dists((samePost>0) & (samePre==0)),hRange);
    difPreHist = histcounts(dists(samePre==0),hRange);
    difPostDifHist = histcounts(dists((samePost==0) & (samePre==0)),hRange);

    samePostHist = histcounts(dists((samePost>0) ),hRange);
    %difPreHist = hist(dists(samePre==0),hRange);
    difPostHist = histcounts(dists((samePost==0) ),hRange);

    
    
%     
%     clf
%     plot(hRange,samePreHist./(difPreHist+samePreHist),'r')
%     hold on
%     plot(hRange,samePostHist./(difPostHist+samePostHist),'b')
%     hold off
%     ylim([0 1])
%     
%    


    
    ratePostHist = samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1));
    ratePostDifHist = samePostDifHist(1:end-1)./(difPostDifHist(1:end-1)+samePostDifHist(1:end-1));
    %plot(hRange(1:end-1),samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1)),'b')
    bar(showRange(1:end-1),ratePostHist,'b')
    hold on
    bar(showRange(1:end-1),ratePostDifHist,'r')
    bar(showRange(1:end-1),[ratePostHist; ratePostDifHist]')

    xlim([-2 100])
    ylim([0 1])
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
    
    
    
%%%%%%%%%%%%%%% Replicate Unbiased %%%%%%
%%    
syn1 = unique(cat(1,synGroup{:}));
syn2 =  unique(cat(1,synGroup{:})); 
    
samePost = zeros(length(bout),length(bout));
samePre = samePost; dists = samePost;useComp = samePost;
for s1 = 1:length(bout)
    for s2 = 1:length(bout)
        
        %samePost(s1,s2) = syn1(s1),1)==edges(syn2(s2),1);
        %samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
        dists(s1,s2) = sqrt((anc(s1,1) - anc(s2,1)).^2 + ...
            (anc(s1,2) - anc(s2,2)).^2 + ...
            (anc(s1,3) - anc(s2,3)).^2);
        useComp(s1,s2) = s1<s2;
        sameComp(s1,s2) = s1==s2;
        
        
        post1 = posts{s1};
        post2 = posts{s2};
        numIntersect = length(intersect(post1,post2));
        atLeast(s1,s2) = numIntersect>0;
        samePost = atLeast;
        
    end
end

symDists = dists;
symDists(sameComp) = inf;
dists(~useComp) = inf;



hRange = [0:6:100];
showRange = (hRange(2:end)+hRange(1:end-1))/2;
%showRange = hRange(1:end-1);

%[samePreHist] = histcounts(dists(samePre>0),hRange);
samePostHist = histcounts(dists(samePost>0),hRange);
%difPreHist = histcounts(dists(samePre==0),hRange);
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
    
    
    
    
    