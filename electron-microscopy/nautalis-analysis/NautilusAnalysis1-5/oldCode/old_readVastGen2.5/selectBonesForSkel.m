function[surfSkel] = selectBonesForSkel(surfSkel);

%%
minLength = 100;


%%
allBones = surfSkel.allBones;
boneLengths = [allBones.length];
longEnough = boneLengths>=minLength;
surfSkel.bones = allBones(longEnough);

% %%
% useBones = zeros(length(allBones),1);
% for i = 1:length(allBones)
%     
%     badBone = allBones(i)
%     
% end
% 

%%

tipList = [allBones.tip];
parentList = [allBones.parent];
checkParents = unique(parentList);
checkParents = checkParents(checkParents>0);
numNodeList = [allBones.numNodes];

for i = 1:length(checkParents)
    tip = checkParents(i);
    tipTarg = find(tipList == tip);
    children = find(parentList==tip);
    nodeNumbers = numNodeList(children);
    numKids(i) = length(children);
    medLength(i) = median(nodeNumbers);
    parentLength(i) = numNodeList(tipTarg);
end



%% Get width


parfor i = 1:length(allBones)
    tip = allBones(i).tip;
    children = find(parentList==tip);
    numKids(i) = length(children);
end




