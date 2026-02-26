function[indi] = checkNodeIndependence(skel,nearestNode,con);



%con = arbor.vox.conMat;
node2branch = skel.node2bones';
vox2branch = node2branch(nearestNode);

selfMat = repmat([1:size(con,1)]', [1 26]);
con2 = con;
con2(con2==0) = selfMat(con2==0);
touchBranchMat = vox2branch(con2); %branch identity of touched voxels


%% find voxels that touch voxels from other branches
touchDif = sum(touchBranchMat ~= repmat(vox2branch,[1 26]),2);

nodeDif = node2branch * 0;
for n = 1:max(nearestNode(:))
    nodeDif(n) = sum(touchDif(nearestNode==n));
end
indiNode = nodeDif==0;
indiBranchNodes = zeros(max(node2branch),1);
for b = 1:max(node2branch)
    indiBranchNodes(b) = sum(indiNode(node2branch==b));
end


if 0
    col = jet(100);
    %col = col(randperm(100),:);
    nodeVal = double(nodeDif>0);
    
    colIDX = ceil(((nodeVal)/max(nodeVal))*99) + 1;
    nodeCol = col(colIDX,:);
    clf
    showBones3D4([],skel,nodeCol);
end


indi.nodes = indiNode;
indi.branches = indiBranchNodes;















