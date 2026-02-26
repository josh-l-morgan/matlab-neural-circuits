function[arbor] = group2nodes(surfVox)

%% get nodes
nodeInd = surfVox.skel.node2surf;
nodeProp = zeros(size(surfVox.subs,1),1);
nodeProp(nodeInd) = 1;

%%

reps = 10;
[surfClose2node] = passMaxAndDist(surfVox,nodeProp,reps);
surfVox.surfClose2node = surfClose2node;

showMaxProp(surfVox.subs,mod(surfClose2node.vox*7777,256));


isSurf = sum(surfVox.conMat>0,2)<20;
nearNum = nodeInd*0;
nearMide = zeros(length(nodeInd),3);
numSurf = nodeInd * 0;
surfRad = nodeInd * 0;
for i = 1:length(nodeInd)
    nearVox = find(surfClose2node.vox == nodeInd(i) );
    nearNum(i) = length(nearVox);
    nearSurf = nearVox(isSurf(nearVox));
    nearSubs = surfVox.subs(nearVox,:);
    nearMid(i,:) = mean(nearSubs,1);
    nearDist = sqrt((nearSubs(:,1) - nearMid(i,1)).^2 + (nearSubs(:,2)-nearMid(i,2)).^2 + ...
        (nearSubs(:,3)-nearMid(i,3)).^2);
    nearMeanDist(i) = mean(nearDist);
    
    %%Cylinder calculation
    if isempty(nearSurf)
        %cylSurf(i) = 1;
        nearSurfSubs = surfVox.subs(nearVox,:);
        surfMid = mean(nearSurfSubs,1);
        nearDist = sqrt((nearSurfSubs(:,1) - surfMid(1)).^2 + (nearSurfSubs(:,2)-surfMid(2)).^2 + ...
            (nearSurfSubs(:,3)-surfMid(3)).^2);
        nearDist = sort(nearDist,'descend');
        
        surfRad(i) = mean(nearDist(1:min(4,length(nearDist)))); %find radius using 4 furthest
             
         mean(nearDist(1:min(4,length(nearDist))))
        numSurf(i) = 0;
    else
    nearSurfSubs = surfVox.subs(nearSurf,:);
    surfMid = mean(nearSurfSubs,1);
    nearDist = sqrt((nearSurfSubs(:,1) - surfMid(1)).^2 + (nearSurfSubs(:,2)-surfMid(2)).^2 + ...
        (nearSurfSubs(:,3)-surfMid(3)).^2);
    surfRad(i) = mean(nearDist);
    numSurf(i) = length(nearSurf);
    
    end

end

%%
cylSurf =  2 * pi * surfRad * surfVox.skel.mergeSize;
cylVol = pi * (surfRad-.5).^2 * surfVox.skel.mergeSize;
cylSurf2vol = cylSurf./(cylVol+cylSurf);
surf2vol = numSurf./nearNum;
    
surfQual = surf2vol-cylSurf2vol;

showMaxProp(surfVox.skel.node2subs,(surfQual+1)*100);
pause(.01)

near.node2surf = nodeInd;
near.num = nearNum;
near.numSurf = numSurf;
near.mid = nearMid;
near.dist = nearMeanDist;
near.surfQual = surfQual;
near.vox = nearVox;


%% decide on bones

for i = 1:length(surfVox.skel.bones)
    boneNodes = surfVox.skel.bones(i).nodes;
    boneQuals = surfQual(boneNodes);
%     showNodes = surfVox.skel.node2bones * 0;
%     allNodes  = showMaxProp(surfVox.skel.node2subs,showNodes+100);
%     showNodes(boneNodes) = 1000;
%     currentNodes  = showMaxProp(surfVox.skel.node2subs,showNodes*1000);
%     subplot(2,1,1)
%     image(allNodes(500:end,700:end) + currentNodes(500:end,700:end))
%     subplot(2,1,2)
%     plot(boneQuals)
%     mean(boneQuals)
    meanBoneSurfQual(i) = mean(boneQuals);
      
    
end


%% Get good bones
boneSurfQualThresh = -.25;
goodBones = meanBoneSurfQual>=boneSurfQualThresh;
skel = surfVox.skel;
clear arbor
arbor.branches = skel.bones(goodBones);
branchNodes = unique([cat(2,arbor.branches.nodes) cat(2,arbor.branches.base) cat(2,arbor.branches.tip) cat(2,arbor.branches.parent)]);
%% make nodes unique vox
node2vox = skel.node2surf(branchNodes);
[uniqueNode2vox ia ic] = unique(node2vox);
uniqueBranchNodes = branchNodes(ia);
skelNode2BranchNode = zeros(max(branchNodes),1);
skelNode2BranchNode(branchNodes) = ic;


%%
node2branch = [];
for i = 1:length(arbor.branches)
    arbor.branches(i).nodes  = skelNode2BranchNode(arbor.branches(i).nodes);
    node2branch = cat(1,node2branch,ones(length(arbor.branches(i).nodes),1)*i);
    arbor.branches(i).tip  = skelNode2BranchNode(arbor.branches(i).tip);
    arbor.branches(i).base  = skelNode2BranchNode(arbor.branches(i).base);
    arbor.branches(i).parent = skelNode2BranchNode(arbor.branches(i).parent);
    arbor.branches(i).edges  = skelNode2BranchNode(arbor.branches(i).edges);
    arbor.branchLengths(i) = sum(arbor.branches(i).edgeLengths);
end    
    
arbor.nodes.node2vox = skel.node2surf(uniqueBranchNodes);
arbor.nodes.node2subs = skel.node2subs(uniqueBranchNodes,:);
arbor.nodes.node2branch = node2branch;
arbor.nodes.seed = skel.surfSeed;
arbor.bridges = skelNode2BranchNode(skel.bridges);

arbor.vox.subs = surfVox.subs;
arbor.vox.conMat = surfVox.conMat;
arbor.vox.conDat = surfVox.conDat;
arbor.vox.seed = surfVox.seedPath.seed;
arbor.vox.segID = surfVox.seedPath.segID;

arbor.offset = skel.offset;
arbor.fsize = skel.fsize;

if length(arbor.nodes.node2vox) ~= length(unique(arbor.nodes.node2vox))
    'NODES ARE NOT UNIQUE'
    pause(2)
end

%% fix bridges
bridges = skel.bridges;
newBridges = bridges * 0;
for i = 1:length(bridges(:));
    bridgePos = skel.node2subs(bridges(i),:);
    dists = sqrt((arbor.nodes.node2subs(:,1)-bridgePos(1)).^2 + ...
        (arbor.nodes.node2subs(:,2)-bridgePos(2)).^2 + (arbor.nodes.node2subs(:,3)-bridgePos(3)).^2);
   newBridges(i) = find(dists == min(dists),1); 
end
arbor.bridges = newBridges;

%% New closest vox

nodeInd = arbor.nodes.node2vox;
nodeProp = zeros(size(arbor.vox.subs,1),1);
nodeProp(nodeInd) = 1;

%% find closest with filtered results

reps = 10;
[surfClose2node] = passMaxAndDist(surfVox,nodeProp,reps);
surfVox.surfClose2node = surfClose2node;

showMaxProp(surfVox.subs,mod(surfClose2node.vox*7777,256));
numVox = nodeInd * 0;
%% Oportunity to grab volume statistics
for i = 1:length(nodeInd) %run each node
    nearVox = find(surfClose2node.vox == nodeInd(i) );
    nearestNode(nearVox) = i;
    numVox(i) = length(nearVox);
end

arbor.vox.nearestNode = nearestNode;
arbor.node.numVox = numVox;















