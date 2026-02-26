function[skel] = simpleBones(surfVox,mergeSize);
%%
if ~exist('mergeSize')
    mergeSize = 5;
end

skel.offset = surfVox.minMax{1};
skel.fsize = surfVox.minMax{2};


subs = surfVox.subs;
for i = 1:3
    subs(:,i) = subs(:,i) - skel.offset(i);
end
inds = sub2ind(skel.fsize,subs(:,1),subs(:,2),subs(:,3));

% 
% inds2surf = sparse(prod(skel.fsize),1);
% inds2surf(inds) = 1:length(inds);

bones = surfVox.surfSkel.bones;

tips = [bones.tip];
bases = [bones.base];
fixedNode = zeros(size(surfVox.surfSkel.prop));
fixedNode(tips) = tips;
fixedNode(bases) = bases;

useBones = find([bones.use]);
for b = 1:length(useBones)
    oldBone = bones(useBones(b);
    newBones(b).tip = oldBone.tip;
    newBones(b).base = oldBone.base;
    newBons(b).parent = oldBone.parent;
    newBons(b).parent = oldBone.tip;
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
       
        
        if makeNode
            countNode = countNode + 1;
            holdSubs = surfVox.subs(holdVox);
            voxDist = sqrt((y-mean(y)).^2 + (x-mean(x)).^2 + (z-mean(z)).^2);
            nearest = find(voxDist == min(voxDist),1);
            pickSub = [y(nearest(1)) x(nearest(1)) z(nearest(1))];
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
    baseSub = indToSub(skel.fsize,oldBone.base);
    edgeLength(countNode) =  sqrt((pickSub(1)-baseSub(1))^2 + ...
        (pickSub(2)-baseSub(2))^2 + (pickSub(3)-baseSub(3))^2 );
    edgeMid(countNode,:) = mean([pickSub;baseSub],1);
    newBones(b).nodes = newNode;
    newBones(b).edgeLengths = edgeLength;
    newBones(b).edgeMids = edgeMid;
    newBones(b).nodeSubs = indToSub(skel.fsize,newNode);
end

skel.bones = newBones;


