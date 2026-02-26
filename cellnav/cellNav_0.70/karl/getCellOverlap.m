%
%Plan
% take a list of N cids. Return an N x N matrix showing the size of overlap
% in DS voxels.
% Input: cidList, dsObj, obi or tis
% steps:
%  for each cell, get a list of object IDs
%  get the voxels for each cid
%  do a convex hull on a zprojection
%  compare via same code as overlap

function outputMat=getCellOverlap(cids,DS,curtis)

%Get all of the voxels for the cells
cidList=cids;
dsObj=DS;
obI=curtis.obI;
allVox={};
allHulls={};
for cidIt=1:length(cidList)
    curCid=obI.cell.name(cidIt);
    obIdList=obI.cell.obIDs{cidIt};
    cidVoxList=[];
    for obIt=1:length(obIdList)
        cidVoxList=[cidVoxList;dsObj(obIdList(obIt)).subs];
        
    end
    allVox{cidIt}=cidVoxList;
    cidZproj=sum(cidVoxList,3);
end

comboList=nchoosek(1:length(cidList),2);
overlapList=zeros(length(comboList),1);
for compIt=1:length(comboList)
    compIt
    compAID=comboList(compIt,1);
    compBID=comboList(compIt,2);
    compAcid=cidList(compAID);
    compBcid=cidList(compBID);
    compAVox=allVox{compAID};
    compBVox=allVox{compBID};
    if length(compAVox)>1&length(compBVox)>1
    overlap=intersect(compAVox,compBVox,'rows');
    overlapList(compIt)=length(overlap(:,1));
    end
end




end