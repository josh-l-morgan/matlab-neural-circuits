function[skel] = simpleBones(surfVox,mergeSize);
%%
clear node2surf node2subs node2bone skel

if ~exist('mergeSize')
    mergeSize = 5;
end

surfSkel = surfVox.surfSkel;
bones = surfSkel.bones;

%skelPred = surfSkel.pred(surfSkel.prop>0);

skel.offset = surfVox.minMax{1};
skel.fsize = surfVox.minMax{2};
skel.surfSeed = surfVox.firstSeed;

nodeID = 0;

%%Make fixed Nodes
tips = [bones.tip];
bases = [bones.base];
fixedNode = zeros(size(surfSkel.prop));
fixedNode(tips) = tips;
fixedNode(bases) = bases;

useBones = find([bones.use]);
for b = 1:length(useBones)
    oldBone = bones(useBones(b));
    
    nodeID = nodeID + 1;
    node2surf(nodeID) = oldBone.tip;
    node2bone(nodeID) = b;
    newBones(b).tip = nodeID;
    
    nodeID = nodeID + 1;
    node2surf(nodeID) = oldBone.base;
    newBones(b).base = nodeID;
    node2bone(nodeID) = b;
    
    
    nodeID = nodeID + 1;
    if oldBone.parent
    node2surf(nodeID) = oldBone.parent;  
    else
            node2surf(nodeID) = oldBone.tip;  
    end
    newBons(b).parent = nodeID;
    node2bone(nodeID) = b;
    
    
    oldNodes = oldBone.nodes;
    
%     nodeIDlist = nodeID+1:nodeID+length(oldNodes);
%     nodeID = max(nodeIDlist);
%     node2surf(nodeIDlist) = oldNodes;
%     node2bone(nodeIDlist) = b;
%     newBones(b).rawNodes = nodeIDlist;
        
    edges = [];
    edgeList = [];
    edgeMid = [];
    newNode = [];
    edgeLength = [];
    countNode = 0;
    
    
    holdVox = [];
    lastNode = [];
    
    numRaw = length(oldNodes);
    fixedNode(oldNodes(end)) = 1;
    %fixedNode = fixedNode * 0;
    for r = 1:numRaw
        tip = oldNodes(r);
        holdVox = [holdVox tip];
        makeNode = 0;
        
        if   fixedNode(tip) % if this vox is fixed
            makeNode = 1;
        elseif fixedNode(min(oldNodes(r+1),numRaw))
             makeNode = 1;
        elseif (length(holdVox)>=mergeSize) % if hold Vox is ful
            makeNode = 1;
        end
        %
        %         elseif surfSkel.pred(tip)<1; % if next node is end of line
        %             makeNode = 1;
        %          elseif fixedNode(surfSkel.pred(tip)); % if next node is fixed
        %             makeNode = 1;
        %         elseif finishedSkel(tip)
        %             makeNode = 1;
        %         end
        
        if makeNode
            countNode = countNode + 1;
            heldSubs = surfVox.subs(holdVox,:);
            heldPos = mean(heldSubs,1);
            voxDist = sqrt((heldSubs(:,1)-heldPos(1)).^2 + (heldSubs(:,2)-heldPos(2)).^2 +...
                (heldSubs(:,3)-heldPos(3)).^2);
            nearest = find(voxDist == min(voxDist),1);
            pickSub = heldSubs(nearest(1),:);
            pickVox = holdVox(nearest(1));
            
            %newNode(countNode) = pickVox;
            nodeID = nodeID + 1;
            newNode(countNode) = nodeID;
            node2surf(nodeID) = pickVox;
            node2bone(nodeID) = b;
            
            
            if countNode>1;
                %                 newEdge = cat(2,lastNode,holdVox(nearest));
                %                 edgeList = cat(1,edgeList,newEdge);
                %                 newPred(lastNode) = holdVox(nearest);
                edgeLength(countNode-1) = sqrt((pickSub(1)-lastSub(1))^2 + ...
                    (pickSub(2)-lastSub(2))^2 + (pickSub(3)-lastSub(3))^2 );
                edgeMid(countNode-1,:) = mean([pickSub;lastSub],1);
                edges(countNode -1,1:2) = [lastNode nodeID];
            end
            
            lastNode = nodeID;
            lastSub = pickSub;
            holdVox = [];
        end %make node
        
        
    end
    baseSub = surfVox.subs(oldBone.base,:);
    edgeLength(countNode) =  sqrt((pickSub(1)-baseSub(1))^2 + ...
        (pickSub(2)-baseSub(2))^2 + (pickSub(3)-baseSub(3))^2 );
    edgeMid(countNode,:) = mean([pickSub;baseSub],1);
    edges(countNode,1:2) = [lastNode newBones(b).base];
    newBones(b).nodes = newNode;
    newBones(b).edges = edges;
    newBones(b).edgeLengths = edgeLength;
    newBones(b).edgeMids = edgeMid;
    %newBones(b).nodeSubs = surfVox.subs(newNode,:);
end


%% assemble skel
skel.bones = newBones;
skel.node2surf = node2surf;
skel.node2bones = node2bone;
skel.node2subs = surfVox.subs(node2surf,:);


