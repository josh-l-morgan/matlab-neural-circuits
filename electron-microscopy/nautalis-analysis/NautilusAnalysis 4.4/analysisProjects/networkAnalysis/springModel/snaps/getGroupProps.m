function[groupPropIDs groupProps] = getGroupProps(groupNodes,checkIDs,checkProp);

if size(checkProp,1) == 1;
    checkProp = checkProp';
end

for i = 1:length(groupNodes)
    gNodes = groupNodes{i};
    [overlap ai bi] = intersect(checkIDs, gNodes);
    groupPropIDs{i} = checkIDs(ai);
    groupProps{i} = checkProp(ai,:);
end
