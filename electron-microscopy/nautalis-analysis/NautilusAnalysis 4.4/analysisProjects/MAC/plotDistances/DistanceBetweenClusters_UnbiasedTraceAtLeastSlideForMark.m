


clear all
%MPN = GetMyDir;
%load('MPN.mat')
MPN = 'C:\Users\jlmorgan\Documents\Data\Segmentations\AndermannChen\export\MAC_merge_mat\';

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

clf

%% bouton perspective

bout = unique(edges(:,2));

clear anc
for b = 1:length(bout)
    targs = edges(:,2)==bout(b);
    bPos = anchors(targs,:);
    anc(b,:) = mean(bPos,1);
    posts{b} = edges(targs,1);
end

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
multMark = {'d','s','*'};
multSize = [55,25,20];
for i = 1:length(bout)
    
    post = posts{i};
    if length(post)==1;
        scatter(anc(i,1),anc(i,2),45,lookupCol(post,:),...
            'filled','markeredgecolor','w')
    else
        for p = 1:length(post)
            scatter(anc(i,1),anc(i,2),multSize(p),lookupCol(post(p),:),...
                multMark{p},'filled','markeredgecolor','w')
        end
    end
    %set(gca,'MarkerEdgeColor',lookupCol(edges(syn1(i))));
    %scatter(anchors(syn1(i),1),anchors(syn1(i),2),lookupCol(edges(syn1(i),1),:))
    %scatter(anchors(syn1(i),1),anchors(syn1(i),2),50,edges(1,1),'o','filled')
    hold on
end
hold off
xlim([140 210])
ylim([220 290])

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

%% show hists
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


%% Sliding Bins

binWidth = 4;
binStep = 0.25;
showRange = [2:binStep:50];

binRad = binWidth/2;
sameDists = dists(samePost>0);
difDists = dists(samePost == 0);

for i = 1:length(showRange);
    
   samePostHist(i) = sum((sameDists>= (showRange(i)-binRad)) & ...
       (sameDists < (showRange(i) + binRad)));
   difPostHist(i) = sum((difDists>= (showRange(i)-binRad)) & ...
       (difDists < (showRange(i) + binRad)));
        
end

subplot(3,1,2)
sumBoth = sum(samePostHist) + sum(difPostHist);

hold off
plot(showRange(1:end-1),samePostHist(1:end-1)/sumBoth,'r','linewidth',3)
hold on
plot(showRange(1:end-1),difPostHist(1:end-1)/sumBoth,'k','linewidth',3)
hold off


xlim([-2 max(showRange)+2])
hold on
subplot(3,1,3)
hold off
%plot(hRange(1:end-1),samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1)),'b')
plot(showRange(1:end-1),samePostHist(1:end-1)./(difPostHist(1:end-1)+samePostHist(1:end-1)),'r','linewidth',3)

xlim([-2 max(showRange)+2])
ylim([0 1])
hold off
ylim([0 1])



return
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
        
        
        
        
        % samePost(s1,s2) = edges(syn1(s1),1)==edges(syn2(s2),1);
        samePre(s1,s2) = edges(syn1(s1),2)==edges(syn2(s2),2);
        dists(s1,s2) = sqrt((anchors(syn1(s1),1) - anchors(syn2(s2),1)).^2 + ...
            (anchors(syn1(s1),2) - anchors(syn2(s2),2)).^2 + ...
            (anchors(syn1(s1),3) - anchors(syn2(s2),3)).^2);
        useComp(s1,s2) = s1<s2;
        sameComp(s1,s2) = s1==s2;
        
        pre1 = edges(syn1(s1),2);
        pre2 = edges(syn2(s2),2);
        post1 = edges(edges(:,2)==pre1,1);
        post2 = edges(edges(:,2)==pre2,1);
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








