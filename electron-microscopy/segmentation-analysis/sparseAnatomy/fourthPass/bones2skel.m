function[skel] = bones2skel(skel);

bones = skel.bones;
skel.fsize = skel.minMax{2};


nodes = [];
edges = [];
edgeMids = [];
edgeLengths = [];
nodeSubs = [];
for i = 1: length(bones)
   
   bNodes =  bones(i).nodes;
   nodes = [nodes bNodes];
   
   bEdges = cat(2,bNodes(1:end-1)',bNodes(2:end)');
   bEdges(end+1,:) = [bNodes(end),bones(i).base];
   edges  = cat(1,edges,bEdges);
   edgeMids = cat(1,edgeMids, bones(i).edgeMids);
   edgeLengths = [edgeLengths bones(i).edgeLengths];
   nodeSubs = cat(1,nodeSubs,bones(i).nodeSubs);
end











skel.tips = [bones.tip];
skel.bases = [bones.base];
skel.node2surf = nodes;
skel.edges = edges;
skel.edgeMids = edgeMids;
skel.edgeLengths = edgeLengths;
skel.nodeSubs = nodeSubs;
skel.surf2node(nodes) = 1:length(nodes);



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