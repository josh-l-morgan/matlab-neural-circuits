function[p] = shortestPathEpos(edges,pos,firstSeed);

%%
rawEdges = edges;
vNum = size(pos,1);
maxEdge = max(edges(:));
edges = sort(edges,2);
edgeInd = sub2ind([maxEdge maxEdge],edges(:,1),edges(:,2));
uEdge = unique(edgeInd);
[y x] = ind2sub([maxEdge maxEdge],uEdge);
edges = [y x];

%% study
tic
if ~exist('firstSeed','var')
    firstSeed = 1;
end

seed = firstSeed;

%% Grab voxels
con2 = edges;


%% find shortest path assuming no loops

dist2seed = zeros(vNum,1);
dist2seed(seed) = 0;
pred = zeros(vNum,1);
segID = zeros(vNum,1);
segID(seed)  = seed;
countBreak = 0;

edgeLength = sqrt((pos(edges(:,1),1)-pos(edges(:,2),1)).^2 + ...
    (pos(edges(:,1),2)-pos(edges(:,2),2)).^2 + ...
    (pos(edges(:,1),3)-pos(edges(:,2),3)).^2);


%%

bigNum = sum(edgeLength)*10000;
dist2seed = dist2seed * 0 + bigNum;
newDist = dist2seed;

dist2seed(firstSeed) = 0'

sumDists = sum(dist2seed);

% hold off
% scatter3(pos(:,1),pos(:,2),pos(:,3),'.','k')
% hold on
% scatter3(pos(firstSeed,1),pos(firstSeed,2),pos(firstSeed,3),'o',...
%     'filled','r')
% pause(.1)
%%
ittNum = 0;
while 1
   
    ittNum = ittNum +1
    
    %%change edge2
    newDist(edges(:,2),1) = dist2seed(edges(:,1),1) + edgeLength;
    betterDist = newDist(edges(:,2))< dist2seed(edges(:,2));
    pred(edges(betterDist,2)) = edges(betterDist,1);
    dist2seed(edges(betterDist,2)) = newDist(edges(betterDist,2));
        
    %change edge1
    newDist(edges(:,1),1) = dist2seed(edges(:,2),1) + edgeLength;
    betterDist = newDist(edges(:,1))< dist2seed(edges(:,1));
    pred(edges(betterDist,1)) = edges(betterDist,2);
    dist2seed(edges(betterDist,1)) = newDist(edges(betterDist,1));
    
    
%     use = dist2seed<bigNum;
%     scatter3(pos(use,1),pos(use,2),pos(use,3),'o','filled','g');
%     pause(.1)
    
    if sumDists == sum(dist2seed)
        break
    end
    sumDists = sum(dist2seed);
end



if 1
scatter3(pos(:,1),pos(:,2),pos(:,3),4,cVal,'filled')
hold on
scatter3(pos(firstSeed,1),pos(firstSeed,2),pos(firstSeed,3),40,'o',...
     'filled','k')
 hold off
 pause(.1)
end

p.pred = pred;
p.dist2seed = dist2seed;
p.seed = firstSeed;




