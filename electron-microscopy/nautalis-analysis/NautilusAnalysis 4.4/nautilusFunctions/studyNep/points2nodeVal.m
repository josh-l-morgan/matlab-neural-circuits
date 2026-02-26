function[rec] = points2nodeVal(pos1,pos2,vals);

if ~exist('vals','var')
    vals = ones(size(pos2,1));
end

posDists = sqrt((pos2(:,1)-pos1(:,1)').^2 + (pos2(:,2)-pos1(:,2)').^2  + (pos2(:,3)-pos1(:,3)').^2);
minDist = min(posDists,[],2);
isNode = minDist * 0;
nodeVal = zeros(size(pos1,1),1);
valNum = zeros(size(pos1,1),1);
for p = 1:size(pos2,1)
    targ = find(posDists(p,:) == minDist(p),1);
    isNode(p) = targ;
    nodeVal(targ) = nodeVal(targ) + vals(p);
    valNum(targ) = valNum(targ) + 1;
    
end

rec.nodeVal = nodeVal;
rec.valNum = valNum;
rec.isNode = isNode;
