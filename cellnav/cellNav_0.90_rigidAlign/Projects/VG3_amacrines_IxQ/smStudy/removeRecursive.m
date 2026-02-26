function[nep] = removeRecursive(nep)

edges = nep.edges;

goodEdge = find(abs((edges(:,1)-edges(:,2)))>0); %find edges with different nodes

nep.edges = nep.edges(goodEdge,:);
if isfield(nep,'edgeRad')
    nep.edgeRad = nep.edgeRad(goodEdge,:);
end
if isfield(nep,'groupEdges')
    nep.groupEdges = nep.groupEdges(goodEdge);
end