function[smallSub] = shrinkSub(sub,downScale);

if ~isempty(sub)
sub = double(sub);
smallSub = ceil(sub/downScale);
maxSub = max(smallSub,[],1)+2;
subInd = sub2ind(maxSub,smallSub(:,1),smallSub(:,2),smallSub(:,3));
uInd = unique(subInd);
[y x z] = ind2sub(maxSub,uInd);
smallSub = cat(2,y,x,z);
else
    smallSub = sub;
end 