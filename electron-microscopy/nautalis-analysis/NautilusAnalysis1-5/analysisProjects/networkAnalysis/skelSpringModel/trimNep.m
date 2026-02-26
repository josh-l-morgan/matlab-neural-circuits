function[nep2] = trimNep(nep,useNodes);


%%
useNodes = intersect(nep.nodes,useNodes);


nep2 = nep;

nodeFields = {'nodes','nodeName','nodePos','nodeType','nodeParent','useNode',...
    'nodeCol','nodePushes','nodeIsPushed','nodeMass','nodeMove','usePos'};
for i = 1:length(nodeFields)
    if isfield(nep,nodeFields{i});
        f = getfield(nep,nodeFields{i});
        if size(f,2)>size(f,1); f = f';end
        nep2 = setfield(nep2,nodeFields{i},f(useNodes,:));
    end
end
nep2.nodes = 1:length(useNodes);



%edges = uniqueEdgesSimple(nep.edges,conNodes); 
[useEdge newEdges] = usedEdges(nep.edges,useNodes); 

edgeFields = {'edges','edgeLengths','edgeType','useEdge','edgeCol',...
    'edgeWeights','allWeights','edgeWidth'};
for i = 1:length(edgeFields)
    if isfield(nep,edgeFields{i});
        f = getfield(nep,edgeFields{i});
        if size(f,2)>size(f,1); f = f';end
        nep2 = setfield(nep2,edgeFields{i},f(useEdge,:));
    end
end

nep2.edges = newEdges;
nep2.nodeFields = nodeFields;
nep2.edgeFields = edgeFields;



