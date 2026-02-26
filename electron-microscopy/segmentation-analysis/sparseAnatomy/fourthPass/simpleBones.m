function[skel] = simpleBones(surfVox,mergeSize);
%%
if ~exist('mergeSize')
    mergeSize = 5;
end

surfSkel = surfVox.surfSkel;
skelPred = surfSkel.pred(surfSkel.prop>0);

bones = surfSkel.bones;


%
% uPred = unique(skelPred);
% noPred = uPred(uPred < 1);
% uPred = uPred(uPred >= 1);
% newPred = zeros(size(surfSkel.prop));

% numPred = hist(skelPred,uPred);
%
%
% fixedNode(surfSkel.tips) = 1;
% fixedNode(uPred(numPred>1)) = 1;
% fixedNod(surfSkel.prop==0) = 1;
% finishedSkel = fixedNode * 0;

tips = [bones.tip];
bases = [bones.base];
fixedNode = zeros(size(surfSkel.prop));
fixedNode(tips) = tips;
fixedNode(bases) = bases;

useBones = find([bones.use]);
for b = 1:length(useBones)
    oldBone = bones(useBones(b));
    newBones(b).tip = oldBone.tip;
    newBones(b).base = oldBone.base;
    newBons(b).parent = oldBone.parent;
    newBones(b).rawNodes = oldBone.nodes;
    oldNodes = oldBone.nodes;
    
    edgeList = [];
    newNode = [];
    edgeLength = [];
    countNode = 0;
    
    holdVox = [];
    lastNode = [];
    
    numRaw = length(oldNodes);
    fixedNode(newBones(b).rawNodes(end)) = 1;
    
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
            newNode(countNode) = pickVox;
            
            if countNode>1;
                %                 newEdge = cat(2,lastNode,holdVox(nearest));
                %                 edgeList = cat(1,edgeList,newEdge);
                %                 newPred(lastNode) = holdVox(nearest);
                edgeLength(countNode-1) = sqrt((pickSub(1)-lastSub(1))^2 + ...
                    (pickSub(2)-lastSub(2))^2 + (pickSub(3)-lastSub(3))^2 );
                edgeMid(countNode-1,:) = mean([pickSub;lastSub],1);
            end
            
            lastNode = pickVox;
            lastSub = pickSub;
            holdVox = [];
        end %make node
        
        
    end
    baseSub = surfVox.subs(oldBone.base,:);
    edgeLength(countNode) =  sqrt((pickSub(1)-baseSub(1))^2 + ...
        (pickSub(2)-baseSub(2))^2 + (pickSub(3)-baseSub(3))^2 );
    edgeMid(countNode,:) = mean([pickSub;baseSub],1);
    newBones(b).nodes = newNode;
    newBones(b).edgeLengths = edgeLength;
    newBones(b).edgeMids = edgeMid;
    newBones(b).nodeSubs = surfVox.subs(newNode,:);
end

% 
% node2surf =  unique(edgeList(:));
% surf2node = zeros(1,size(surfVox.subs,1));
% surf2node(node2surf) = 1:length(node2surf);
% 
% 
% skel.nodes = 1:length(node2surf);
% skel.edges = surf2node(edgeList);
% skel.tips = surf2node(surfSkel.tips);
% skel.node2surf = node2surf;
% skel.surf2node = surf2node;
% skel.nodeSubs = surfVox.subs(node2surf,:);
% skel.prop = zeros(1,size(surfVox.subs,1));
% skel.prop(node2surf) = 1;
% skel.nodeCon = hist(skel.edges(:),skel.nodes);
% 
% nodePred = newPred(node2surf);
% nodePred(nodePred>0) = surf2node(nodePred(nodePred>0));
% skel.pred = nodePred;

skel.bones = newBones;


