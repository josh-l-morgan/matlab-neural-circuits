function [predList,usedEdgeList]=getBranch(allEdges,curNode,usedEdgeList)
    prevEdge=0;
    nextEdge=find(any(allEdges==curNode,2));
    nextEdge=nextEdge(nextEdge~=prevEdge);
    if length(nextEdge)==1
        edgeNodes=allEdges(nextEdge,:);
        nextNode=edgeNodes(edgeNodes~=curNode);
    else
        recurse
    end


    uniqueNodes=curSM.nep.nodes;
    nodeCounts=tabulate(allEdges(:));
    tipIDs=nodeCounts(nodeCounts(:,2)==1);
    forkIDs=nodeCounts(nodeCounts(:,2)>2);



end