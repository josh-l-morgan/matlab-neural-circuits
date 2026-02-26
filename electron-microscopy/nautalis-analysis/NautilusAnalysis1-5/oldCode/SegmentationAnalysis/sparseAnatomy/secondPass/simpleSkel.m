function[skel] = simpleSkel(surfVox,surfSkel,mergeSize);
%%
if ~exist('mergeSize')
    mergeSize = 5;
end

skelPred = surfSkel.pred(surfSkel.prop>0);
uPred = unique(skelPred);
numPred = hist(skelPred,uPred);

fixedNode = zeros(size(surfSkel.prop));
fixedNode(surfSkel.tips) = 1;
fixedNode(uPred(numPred>1)) = 1;

finishedSkel = fixedNode * 0;


edgeList = [];
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
        elseif fixedNode(surfSkel.pred(tip)); % if next node is fixed
            makeNode = 1;
        elseif surfSkel.pred(tip)<1; % if next node is end of line
            makeNode = 1;
        elseif finishedSkel(tip)
            makeNode = 1;
        end
        
        if makeNode
            
            heldSubs = surfVox.subs(holdVox,:);
            heldPos = mean(heldSubs,1);
            voxDist = sqrt((heldSubs(:,1)-heldPos(1)).^2 + (heldSubs(:,2)-heldPos(2)).^2 +...
                (heldSubs(:,3)-heldPos(3)).^2);
            nearest = find(voxDist == min(voxDist),1);
            pickSub = heldSubs(nearest(1),:);
            if ~isempty(lastNode)
                newEdge = cat(2,lastNode,holdVox(nearest));
                edgeList = cat(1,edgeList,newEdge);
            end
            
            lastNode = holdVox(nearest);
            holdVox = [];
        end
        
        if finishedSkel(tip), break, end % stop when you reach
        
        tip = surfSkel.pred(tip);
    end
end

skel.edges = edgeList;
skel.tips = surfSkel.tips;
skel.nodes = unique(edgeList(:));
skel.prop = zeros(1,size(surfVox.subs,1));
skel.prop(skel.nodes) = 1;
skel.nodeCon = hist(edgeList(:,2),skel.nodes);

%%Nodes are popping up too many times in edge list.



