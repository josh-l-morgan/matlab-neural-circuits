 function cidList=type2cid(typ,subTyp,curtis)
 cidList={};
for i=1:length(typ)
    curTypStr=typ{i};
    curSubTypStr=subTyp{i};
    curTypNum=find(ismember(curtis.cells.type.typeNames,curTypStr));    
    typIDs=find(curtis.cells.type.typeID==curTypNum);
    if string(curSubTypStr)=="all"
        subTypIDs=1:length(curtis.cells.type.subTypeID);
    else
        curSubTypNum=find(ismember(curtis.cells.type.subTypeNames{curTypNum},curSubTypStr));
        subTypIDs=find(curtis.cells.type.subTypeID==curSubTypNum);
    end
    cidIDs=intersect(typIDs,subTypIDs);
    cidList{i}=curtis.cells.cids(cidIDs)';
end
end