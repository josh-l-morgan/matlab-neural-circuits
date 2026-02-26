function[p] = shortestPathToSeed(edges,edgeLength,firstSeed,pos);


if exist('pos','var')
    show = 1;
else
    show = 0;
end

seed = firstSeed;
%%
rawEdges = edges;
vNum = max(edges(:));
maxEdge = max(edges(:));
edges = sort(edges,2);
edgeInd = sub2ind([maxEdge maxEdge],edges(:,1),edges(:,2));
uEdge = unique(edgeInd);
[y x] = ind2sub([maxEdge maxEdge],uEdge);
edges = [y x];


%% find shortest path assuming no loops

dist2seed = zeros(vNum,1);
dist2seed(seed) = 0;
pred = zeros(vNum,1);

bigNum = sum(edgeLength)*10000;
dist2seed = dist2seed * 0 + bigNum;
newDist = dist2seed;

dist2seed(firstSeed) = 0;

sumDists = sum(dist2seed);

ittNum = 0;


if show
    hold off
    scatter3(pos(:,1),pos(:,2),pos(:,3),'.','k')
    hold on
    scatter3(pos(firstSeed,1),pos(firstSeed,2),pos(firstSeed,3),'o',...
        'filled','r')
    pause(.1)
end
%


while ittNum < (vNum*2)
    
    ittNum = ittNum +1;
    
    %%change edge2
    newDist(edges(:,2),1) = dist2seed(edges(:,1),1) + edgeLength;
    betterDist = newDist(edges(:,2))< dist2seed(edges(:,2));
    pred(edges(betterDist,2)) = edges(betterDist,1);
    dist2seed(edges(betterDist,2)) = newDist(edges(betterDist,2));
    sum(betterDist) 
    
    %change edge1
    newDist(edges(:,1),1) = dist2seed(edges(:,2),1) + edgeLength;
    betterDist = newDist(edges(:,1))< dist2seed(edges(:,1));
    pred(edges(betterDist,1)) = edges(betterDist,2);
    dist2seed(edges(betterDist,1)) = newDist(edges(betterDist,1));
    sum(betterDist)
    
    %%test
    
    isClose = find(dist2seed<bigNum);
    hits = []; hitFrom = [];
    for i = 1:length(isClose)
        isFrom = edges(edges(:,1) == isClose(i),2);
        isFrom2 = edges(edges(:,2) == isClose(i),1);
        isCon = [isFrom; isFrom2]
        hits = [hits; isCon];
        hitFrom = [hitFrom; isCon*0+isClose(i)];
    end
    [hitFrom hits]
    setdiff(hits,isClose)

    
    
    if show
        use = dist2seed<bigNum;
        scatter3(pos(use,1),pos(use,2),pos(use,3),'o','filled','g');
        pause(.1)
    end
    
    
    
    if sumDists == sum(dist2seed)
        'no more changes'
        ittNum
        break
    end
    sumDists = sum(dist2seed);
end




if show
    scatter3(pos(:,1),pos(:,2),pos(:,3),4,cVal,'filled')
    hold on
    scatter3(pos(firstSeed,1),pos(firstSeed,2),pos(firstSeed,3),40,'o',...
        'filled','k')
    hold off
    pause(.1)
end





dist2seed(dist2seed == bigNum) = inf;
p.pred = pred;
p.dist2seed = dist2seed;
p.seed = firstSeed;




