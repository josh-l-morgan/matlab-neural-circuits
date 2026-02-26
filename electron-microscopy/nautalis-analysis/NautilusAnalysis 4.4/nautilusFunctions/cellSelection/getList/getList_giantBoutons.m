function[isBig] = getList_giantBoutons(MPN)

%MPN = GetMyDir
load([MPN 'obI.mat'])
% 
% allEdges = obI.nameProps.edges(:,[2 1]);
% targCells = [108 201 109 903 907];
% conTo = makeConTo(obI,targCells);

inSpine = obI.nameProps.inSpine;

preID = [inSpine.pre];
postID = [inSpine.post];
inNum = [inSpine.in];

minCon = 1;
preNames = unique(preID);
histPre = hist(preID,preNames)
preNames = preNames(histPre>=minCon);
preNames = preNames(preNames>0);
postNames = unique(postID);
histPost = hist(postID,postNames);
postNames = postNames(histPost>=minCon);
postNames = postNames(postNames>0);



for i = 1:length(preNames)
    foundC = ((preID == preNames(i))  & (inNum >= 2));
    preCount(i) =sum(foundC);
end


for p = 1:length(postNames)
    foundC = ( (postID == postNames(p)) & (inNum >=2));
    postCount(p) =sum(foundC);
end

isBig = [preNames(preCount>0) postNames(postCount>0)];



