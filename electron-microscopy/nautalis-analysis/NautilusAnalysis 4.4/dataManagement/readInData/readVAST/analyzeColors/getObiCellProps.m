function[obI] = getObiCellProps(obI)


legacy = 1;

cellIDs = obI.nameProps.cellNum;
checkIDs = unique(cellIDs(cellIDs>0));

if isempty(checkIDs) & legacy
    cellIDs = obI.nameProps.tag.firstNum;
    checkIDs = unique(cellIDs(cellIDs>0));
end
    
cellCount = 0;
for i = 1:length(checkIDs)
    obIDs = find((cellIDs == checkIDs(i)) & obI.nameProps.ofID);
    mainID = intersect(obIDs,find(obI.nameProps.tag.cell));
    anyIDs =  find((cellIDs == checkIDs(i)) );
    
    if length(mainID)>1
        obI.nameProps.names(mainID)
        mainID
        'too many main IDs'
        mainID = mainID(1);
    elseif isempty(mainID)
        if ~isempty(obIDs)
            mainID = obIDs(1);
        else
            mainID = anyIDs(1);
        end
    end
    
    if ~isempty(mainID)
        cellCount = cellCount+1;
        obI.cell.obIDs{cellCount} = obIDs;
        obI.cell.name(cellCount) = checkIDs(i);
        obI.cell.mainID(cellCount) = mainID;
        obI.cell.label{cellCount} = obI.nameProps.names{mainID};
        obI.cell.anchors(cellCount,:) =  obI.colStruc.anchors(mainID,:);
    end
    
    
    %%get seed
    cellSeedOb  =  intersect(obI.cell.obIDs{cellCount}, obI.nameProps.tag.seed);
    cellCBOb =  intersect(obI.cell.obIDs{cellCount}, obI.nameProps.tag.soma);
    obI.cell.seedID{cellCount} = cellSeedOb;
    ob.cell.cbID{cellCount} = cellCBOb;
    
    
    if ~isempty(cellSeedOb)
        obI.cell.anchors(cellCount,:) =  obI.colStruc.anchors(cellSeedOb(1),:);
    elseif ~isempty(cellCBOb)
        obI.cell.anchors(cellCount,:) =  obI.colStruc.anchors(cellCBOb(1),:);
    end
    
    
end










