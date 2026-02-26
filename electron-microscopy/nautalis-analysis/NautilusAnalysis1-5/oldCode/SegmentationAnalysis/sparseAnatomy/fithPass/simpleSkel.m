function[skel] = simpleSkel(surfVox,mergeSize);
%%
if ~exist('mergeSize')
    mergeSize = 5;
end

surfSkel = surfVox.surfSkel;
skelPred = surfSkel.pred(surfSkel.prop>0);
uPred = unique(skelPred);
noPred = uPred(uPred < 1);
uPred = uPred(uPred >= 1);
newPred = zeros(size(surfSkel.prop));

numPred = hist(skelPred,uPred);

fixedNode = zeros(size(surfSkel.prop));
fixedNode(surfSkel.tips) = 1;
fixedNode(uPred(numPred>1)) = 1;
fixedNod(surfSkel.prop==0) = 1;
finishedSkel = fixedNode * 0;


edgeList = [];
countNode = 0;
for t = 1:length(surfSkel.tips)
    
    tip = surfSkel.tips(t);
    holdVox = [];
    lastNode = [];
    
    for r = 1:length(finishedSkel)
        if tip<1,break,end
            holdVox = [holdVox tip];
        makeNode = 0;
        if   fixedNode(tip) % if this vox is fixed
            makeNode = 1;
        elseif (length(holdVox)>=mergeSize) % if hold Vox is ful
            makeNode = 1;
        elseif surfSkel.pred(tip)<1; % if next node is end of line
            makeNode = 1;
         elseif fixedNode(surfSkel.pred(tip)); % if next node is fixed
            makeNode = 1;           
        elseif finishedSkel(tip)
            makeNode = 1;
        end
        
        if makeNode
            countNode = countNode + 1;
            heldSubs = surfVox.subs(holdVox,:);
            heldPos = mean(heldSubs,1);
            voxDist = sqrt((heldSubs(:,1)-heldPos(1)).^2 + (heldSubs(:,2)-heldPos(2)).^2 +...
                (heldSubs(:,3)-heldPos(3)).^2);
            nearest = find(voxDist == min(voxDist),1);
            pickSub = heldSubs(nearest(1),:);
            if ~isempty(lastNode)
                newEdge = cat(2,lastNode,holdVox(nearest));
                edgeList = cat(1,edgeList,newEdge);
                newPred(lastNode) = holdVox(nearest);
            end
            
            lastNode = holdVox(nearest);
            holdVox = [];
        end
        
        if finishedSkel(tip), break, end % stop when you reach
        finishedSkel(tip) = 1;
        tip = surfSkel.pred(tip);
    end
end


node2surf =  unique(edgeList(:));
surf2node = zeros(1,size(surfVox.subs,1));
surf2node(node2surf) = 1:length(node2surf);


skel.nodes = 1:length(node2surf);
skel.edges = surf2node(edgeList);
skel.tips = surf2node(surfSkel.tips);
skel.node2surf = node2surf;
skel.surf2node = surf2node;
skel.nodeSubs = surfVox.subs(node2surf,:);
skel.prop = zeros(1,size(surfVox.subs,1));
skel.prop(node2surf) = 1;
skel.nodeCon = hist(skel.edges(:),skel.nodes);

nodePred = newPred(node2surf);
nodePred(nodePred>0) = surf2node(nodePred(nodePred>0));
skel.pred = nodePred;


