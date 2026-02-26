function[sNode] = paired(nPairs,node)
%%find all nodes connected to node by a paired list

% 

sNode = node;

for i = 1:size(nPairs,1)
    startLength = length(sNode);
    tNode = [];
    for s = 1:length(sNode)
        [p n] = find(nPairs == sNode(s));
        t = ~(n-1)+1;
        tNode = [tNode; nPairs(sub2ind(size(nPairs),p,t))];
    end
    sNode = [sNode;tNode];
    sNode = unique(sNode);
    stopLength = length(sNode);
    if startLength == stopLength, break, end

end