%calculateOverlap

load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\dsObj.mat');
load('Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Merge\obI.mat');

cidList=obI.cell.name;
allVox={};
for cidIt=1:length(cidList)
    curCid=obI.cell.name(cidIt);
    obIdList=obI.cell.obIDs{cidIt};
    cidVoxList=[];
    for obIt=1:length(obIdList)
        cidVoxList=[cidVoxList;dsObj(obIdList(obIt)).subs];
        
    end
    allVox{cidIt}=cidVoxList;
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