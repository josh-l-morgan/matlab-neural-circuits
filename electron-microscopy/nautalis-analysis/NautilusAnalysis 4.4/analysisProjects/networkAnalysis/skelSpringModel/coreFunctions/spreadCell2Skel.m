function[p] = spreadCell2Skel(nep,prop,aveReach);

if ~exist('prop','var')
    prop =  nep.nodeCol;
end
if ~exist('aveReach','var')
    aveReach = 100;
end

useProp = prop(:,1) * 0;
pnum = size(prop,2); %number of properties
edges = nep.edges;

%% Shift cellProp to skel nodes. 
nodeType = nep.nodeType;
cellNodes = find(nodeType==1);
synNodes = find(nodeType == 2);
skelNodes = find(nodeType == 3);

%%transmit cell prop to skeleton
for c = 1:length(cellNodes)
    cellProp = prop(cellNodes(c),:);
    synTarg = edges(edges(:,2) == cellNodes(c),1);
    prop(synTarg,:) = repmat(cellProp,[length(synTarg) 1]);
    for s = 1:length(synTarg)
        skelTarg = edges(edges(:,2) == synTarg(s),1);
        prop(skelTarg,:) = repmat(cellProp,[length(skelTarg) 1]);
        useProp(skelTarg) = 1;
    end
end


%% average skel col
isskel = nep.edgeType == 3;

for i = 1:length(skelNodes)
%     scatter3(nep.nodePos(skelNodes,1),nep.nodePos(skelNodes,2),nep.nodePos(skelNodes,3),'.','k')
%     hold on
   
    start = skelNodes(i);
    
    lastNode = start;
    lastLength = 0;
    useEdge = ones(size(edges,1),1);
    flipEdge = [useEdge*2 useEdge];
    storeNodes = []; storeEdges = []; storeLengths = [];
    while 1
    %for spr = 1:5; 
    nextNodes = [];
        nextLength = [];
        nextEdges = [];
        newNum = 0;
        for n  = 1:length(lastNode)
           hitEdgeN = ismember(edges,lastNode) & [useEdge useEdge] & [isskel isskel];
           hitEdges = sum(hitEdgeN,2)>0;    
           hitNodes = edges([hitEdgeN(:,2) hitEdgeN(:,1)]);
           numHit = length(hitNodes);
           nextLength(newNum+1:newNum +numHit)  = nep.edgeLengths(hitEdges)+lastLength(n);
           nextNodes(newNum+1:newNum +numHit) = hitNodes;
           nextEdges(newNum+1:newNum +numHit) = find(hitEdges);
           newNum = newNum+numHit;
           useEdge(hitEdges) = 0;
          
        end
        
        lastLength = nextLength(nextLength<=aveReach);
        lastNode = nextNodes(nextLength<=aveReach);
        
%         scatter3(nep.nodePos(lastNode,1),nep.nodePos(lastNode,2),nep.nodePos(lastNode,3),'o','filled','b')
% pause
        
        storeEdges = [storeEdges nextEdges(nextLength<=aveReach)];
        
        storeNodes = [storeNodes lastNode];
        storeLengths = [storeLengths nextLength(nextLength<=aveReach)];
        
        if isempty(lastNode)
            break
        end
        
    end
    
    
    aveNodes{i} = storeNodes;
    usedEdge{i} = unique(storeEdges);
    nodeDists{i} = storeLengths;
    
%     scatter(nep.nodePos(storeNodes,1),nep.nodePos(storeNodes,2),'o','filled','y')
% hold off
% pause(.01)
    
end

%% Color skeleton

newProp = zeros(length(prop),1);
p.mean = newProp; p.max = newProp; p.sum = newProp; p.length = newProp;
p.aveNodes = aveNodes;
p.usedEdge = usedEdge;
p.nodeDists = nodeDists;

for i = 1:length(skelNodes)
    
    [sampNodes ia ib] = intersect(find(useProp),aveNodes{i});
    sampDists = nodeDists{i}(ib);
    %sampNodes = aveNodes{i};
    if isempty(sampNodes)
        newProp(skelNodes(i),:) = zeros(1,pnum);
    else
        getProp = prop(sampNodes,:);
        
        p.mean(skelNodes(i),:) = mean(getProp,1);
        p.max(skelNodes(i),:) = max(getProp);
        p.sum(skelNodes(i),:) = sum(getProp);
        p.hits(skelNodes(i),:) = length(getProp);
        p.length(skelNodes(i),:) =  sum(nep.edgeLengths(usedEdge{i}));
        p.foundNodes{skelNodes(i)} = sampNodes;
        p.foundDists{skelNodes(i)} = sampDists;
        p.minDist(skelNodes(i)) = min(sampDists);
    end
        
end



%  scatter(nep.nodePos(skelNodes,1),nep.nodePos(skelNodes,2),'.','k')
% hold on
% for i = 1:length(skelNodes)
%      scatter(nep.nodePos(skelNodes(i),1),nep.nodePos(skelNodes(i),2),'o',...
%          'markerfacecolor',newProp(skelNodes(i),:),'markeredgecolor','none');
%      pause(.001)
% end
% hold off



%%







