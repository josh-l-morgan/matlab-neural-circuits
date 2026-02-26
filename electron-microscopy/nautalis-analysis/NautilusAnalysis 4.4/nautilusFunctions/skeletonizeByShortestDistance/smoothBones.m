function[skel] = smoothBones(surfVox,mergeSize);
%%
clear node2surf node2subs node2bone skel nodePos node2surf node2bone newBones

if ~exist('mergeSize')
    mergeSize = 5;
end


surfSkel = surfVox.surfSkel;
bones = surfSkel.bones;

%skelPred = surfSkel.pred(surfSkel.owner>0);

skel.offset = surfVox.minMax{1};
skel.fsize = surfVox.minMax{2};
skel.surfSeed = surfVox.firstSeed;
skel.mergeSize = mergeSize;

nodeID = 0;

%%Make fixed Nodes
tips = [bones.tip];
bases = [bones.base];
fixedNode = zeros(size(surfSkel.owner));
%fixedNode(tips) = tips;
fixedNode(bases) = bases;
%fixedNode = fixedNode+1;
useBones = find([bones.use]);
for b = 1:length(useBones)
    oldBone = bones(useBones(b));
    
    %%Make tip
    nodeID = nodeID + 1;
    node2surf(nodeID) = oldBone.tip;
    node2bone(nodeID) = b;
    newBones(b).tip = nodeID;
    nodePos(nodeID,:) = surfVox.subs(oldBone.tip,:);
    
    %%Make base
    nodeID = nodeID + 1;
    node2surf(nodeID) = oldBone.base;
    newBones(b).base = nodeID;
    node2bone(nodeID) = b;
    nodePos(nodeID,:) = surfVox.subs(oldBone.base,:);
    
    %%Use parent????????????
    nodeID = nodeID + 1;
    if oldBone.parent
        node2surf(nodeID) = oldBone.parent;
    else
        node2surf(nodeID) = oldBone.tip;
    end
    newBones(b).parent = nodeID;
    node2bone(nodeID) = b;
    nodePos(nodeID,:) = surfVox.subs(node2surf(nodeID),:);

    
    oldNodes = [oldBone.tip oldBone.nodes oldBone.base];
    
    newNode = [];
    countNode = 0;
    holdVox = [];
    lastNode = [];
    
    
    
    numRaw = length(oldNodes);
    if numRaw>2
        for r = 2:numRaw-1
            countNode = countNode + 1;
            nodeID = nodeID + 1;
            %holdVox = oldNodes(max(1,r-mergeSize):min(r+mergeSize,length(oldNodes)));
            isfixed = 1; %dont move node by default
            if ((r-mergeSize)>0) & ((r+ mergeSize)<length(oldNodes));
              holdVox = oldNodes(r-mergeSize:r+mergeSize);
              isfixed = sum(fixedNode(holdVox));
            else
              holdVox = oldNodes(r);
            end
            if isfixed %use original position if spread includes a fixed (base) node
                holdVox = oldNodes(r);
            end
            heldSubs = surfVox.subs(holdVox,:);
            heldPos = mean(heldSubs,1);
            voxDist = sqrt((heldSubs(:,1)-heldPos(1)).^2 + (heldSubs(:,2)-heldPos(2)).^2 +...
                (heldSubs(:,3)-heldPos(3)).^2);
            nearest = find(voxDist == min(voxDist),1);
            pickSub = heldSubs(nearest(1),:);
            pickVox = holdVox(nearest(1));
            newNode(countNode) = nodeID;
            node2surf(nodeID) = pickVox;
            node2bone(nodeID) = b;
            nodePos(nodeID,:) = heldPos;
            
        end %make node
        
        newRun = [newBones(b).tip newNode newBones(b).base];
    else
        newRun = [newBones(b).tip newBones(b).base];
    end
    
    
    edges = [newRun(1:end-1)' newRun(2:end)'];
    edgeLength = sqrt((nodePos(edges(:,1),1) - nodePos(edges(:,2),1)).^2 + ...
        (nodePos(edges(:,1),2) - nodePos(edges(:,2),2)).^2 + ...
        (nodePos(edges(:,1),3) - nodePos(edges(:,2),3)).^2);
    edgeMid = (nodePos(edges(:,1),:) + nodePos(edges(:,2),:))/2;
    
    
    newBones(b).nodes = newNode;
    newBones(b).edges = edges;
    newBones(b).edgeLengths = edgeLength;
    newBones(b).edgeMids = edgeMid;
    newBones(b).nodePos = nodePos(newNode,:);
    %newBones(b).nodeSubs = surfVox.subs(newNode,:);
    
    if length(newBones(b).nodes) ~= length(oldBone.nodes)
       'bark'
       
    end
end


%% assemble skel
skel.bones = newBones;
skel.node2surf = node2surf;
skel.node2bones = node2bone;
skel.node2subs = surfVox.subs(node2surf,:);
skel.nodePos = nodePos;


