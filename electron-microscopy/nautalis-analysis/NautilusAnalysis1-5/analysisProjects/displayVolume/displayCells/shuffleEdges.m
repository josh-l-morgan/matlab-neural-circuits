function[newEdge] = shuffleEdges(edges)

%%randomly shuffle position of one side of an edge list.

shiftPost = randperm(size(edges,1));
newEdge(:,1) = edges(shiftPost,1);
newEdge(:,2) = edges(:,2);