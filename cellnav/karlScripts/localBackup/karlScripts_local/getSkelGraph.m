function [predList,usedNodeList]=getSkelGraph(allEdges,curNode,usedNodeList,predList)
%predList=[curNode 0];
while true
    nextEdge=find(any(allEdges==curNode,2));
    edgeNodes=allEdges(nextEdge,:);
    nextNode=edgeNodes(edgeNodes~=curNode&~ismember(edgeNodes,usedNodeList));
    usedNodeList=[usedNodeList;curNode];
    if isempty(nextNode)
        break
    end
    for j=1:length(nextNode)
        predList=[predList;[nextNode(j) curNode]];
    end
    if length(nextNode)>1
        for k=1:length(nextNode)
            [predList,usedNodeList]=getSkelGraph(allEdges,nextNode(k),usedNodeList,predList);
        end
    else
        curNode=nextNode;
    end
end
end