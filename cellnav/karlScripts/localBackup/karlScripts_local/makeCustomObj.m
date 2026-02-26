function objList=makeCustomObj(includeStrings,excludeStrings,tis)
% This will make you a custom object (maybe also FV) from segmentation
% objects with the inclusion and exclusion strings specified.
% Example:
% makeCustomObj({["1216","syn"],["1197","syn"]},{["probably"]},tis)
% ( ( 1216 & syn ) | ( 1197 & syn ) ) & ~probably
% includes '1216' AND 'syn', or '1197' AND 'syn', but nothing with 'probably'
% good luck
allNams=tis.obI.nameProps.names';
allIncList={};
for curIncIt=1:length(includeStrings)
    curIncCell=includeStrings{curIncIt};
    curList={};
    for curStrIt=1:length(curIncCell)
        curStr=curIncCell(curStrIt);
        curList{curStrIt}=find(contains(allNams,curStr));
    end
    finList=curList{1};
    for curIt=2:length(curList)
        finList=intersect(finList,curList{curIt});
    end
    finalIncLists{curIncIt}=finList;
end
allIncList=vertcat(finalIncLists{:});
finalExcLists={};
for curExcIt=1:length(excludeStrings)
    curIncCell=excludeStrings{curExcIt};
    curList={};
    for curStrIt=1:length(curIncCell)
        curStr=curIncCell(curStrIt);
        curList{curStrIt}=find(contains(allNams,curStr));
    end
    finList=curList{1};
    if length(curList)>1
        for curIt=2:length(curList)
            finList=intersect(finList,curList{curIt});
        end
    end
    finalExcLists{curExcIt}=finList;
end
allExcList=vertcat(finalExcLists{:});
objList=setdiff(allIncList,allExcList);



