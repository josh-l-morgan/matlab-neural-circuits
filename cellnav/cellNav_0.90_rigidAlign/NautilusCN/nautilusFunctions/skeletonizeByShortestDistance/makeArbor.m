function[arbor] = makeArbor(surfVox)

clear arbor
skel = surfVox.skel;

arbor.skel = skel;
arbor.subs = surfVox.subs - repmat(skel.offset,[size(surfVox.subs,1) 1])+1;

arbor.offset = skel.offset;
arbor.vox.subs = surfVox.subs;
arbor.vox.conMat = surfVox.conMat;
arbor.vox.conDat = surfVox.conDat;
arbor.vox.voxBridges = surfVox.voxBridges;
arbor.vox.fv = surfVox.voxFV;


%%Check Redundant
bPos = round(skel.nodePos * 1000);
maxPos = max(bPos,[],1);
oldInd = sub2ind(maxPos,bPos(:,1),bPos(:,2),bPos(:,3));
[uInd ia ic] = unique(oldInd);

%%new Nodes
arbor.nodeNum = length(uInd);
arbor.nodes.nodes = 1:arbor.nodeNum;
arbor.nodes.pos = skel.nodePos(ia,:);
arbor.nodes.rad = skel.nodeRad(ia);
arbor.nodes.node2surf = skel.node2surf(ia);

%%New edges
edges = cat(1,skel.bones.edges,skel.bridges);
edges = ic(edges);
arbor.edges = edges;
arbor.edgeProps.rad = mean(arbor.nodes.rad(edges),2);

s1 = arbor.nodes.pos(edges(:,1),:);
s2 = arbor.nodes.pos(edges(:,2),:);
arbor.edgeProps.pos = (s1 + s2)/2;
arbor.edgeProps.length = sqrt((s1(:,1)-s2(:,1)).^2 + ...
    (s1(:,2)-s2(:,2)).^2 + (s1(:,3)-s2(:,3)).^2);



%% legacy fields
arbor.branches = skel.bones;











