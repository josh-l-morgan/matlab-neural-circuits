function[tipShape] = getTipShape(surfVox,voxEdges);

minRegSize = 100;

regIDs = unique(voxEdges(:,1));
regSize = hist(voxEdges(:,1),regIDs);
bigReg = find(regSize>minRegSize);

onEdge = voxEdges(:,1) ~= voxEdges(:,2);
numBig = length(bigReg);



tipLength = zeros(1,size(voxEdges,1));
edge2cent = tipLength;

for i = 1:numBig
    %disp(sprintf('%d of %d',i,numBig))
    ID = regIDs(bigReg(i));
    idPos = voxEdges(:,1) == ID;
    iSubs = surfVox.subs(idPos,:);
    iNum = size(iSubs,1);
    iEdge = onEdge(idPos);
    edgeSubs = iSubs(iEdge,:);
    edgeNum = sum(iEdge);
    %notEdge(bigReg(i)) = (iNum-edgeNum)/iNum;
    
    midAll = mean(iSubs,1);
    midEdge = mean(edgeSubs,1);
    tip = surfVox.subs(ID,:);
    tipLength(ID) = sqrt((tip(1) - midEdge(1))^2 + (tip(2) - midEdge(2))^2 + ...
        (tip(3) - midEdge(3))^2);
    edge2cent(ID) = sqrt((midAll(1) - midEdge(1))^2 + (midAll(2) - midEdge(2))^2 + ...
        (midAll(3) - midEdge(3))^2);
    
end

tipShape.length = tipLength;
tipShape.edge2cent = edge2cent;

