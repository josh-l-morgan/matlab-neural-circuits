function[con allPreCells allPostCells] = edge2con(edges)


[sortPre preIdx] = sort(edges(:,2));
[sortPost postIdx] = sort(edges(:,1));

allPreCells = unique(edges(:,2));
allPostCells = unique(edges(:,1));

lookupPre(allPreCells + 1) = 1:length(allPreCells);
lookupPost(allPostCells + 1) = 1:length(allPostCells);

newPostEdge = lookupPost(edges(:,1)+1);
newPreEdge = lookupPre(edges(:,2)+1);

conSize = [length(allPreCells) length(allPostCells)];
con = zeros(conSize);

edgeInd = sub2ind(conSize,newPreEdge,newPostEdge);
uInd = unique(edgeInd);
hInd = hist(edgeInd,uInd);
con(uInd) = hInd;

