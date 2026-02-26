function[edges pos] = doubleEdge(edges,pos,rep)

%%doubles sampling of edges without changing possitions

if ~exist('rep','var')
    rep = 1;
end

for i = 1:rep
    newPos = mean(cat(3,pos(edges(:,1),:),pos(edges(:,2),:)),3);
    newPosNum = [size(pos,1) + 1 : size(pos,1) + size(newPos,1)]';
    edges = cat(1,[edges(:,1) newPosNum],[newPosNum edges(:,2)]);
    pos = cat(1,pos,newPos);
end