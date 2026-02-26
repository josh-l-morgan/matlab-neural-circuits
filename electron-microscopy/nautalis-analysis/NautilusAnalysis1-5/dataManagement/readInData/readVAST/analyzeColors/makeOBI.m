function[obI] = makeOBI(TPN,MPN)

%MPN =       'D:\LGNs1\Segmentation\VAST\S8\joshm\export+14+04+27_mat\'

%fileName =  [MPN 'S8_joshmHome_14+05+28_export.txt'];
if ~exist('fileName','var')
textDir = dir([TPN '*.txt']);
    TFN = textDir(1).name;
    
    %[TFN TPN] = uigetfile('.txt');
fileName = [TPN TFN];
end


obI.colStruc = readVastColors(fileName);
obI.nameProps = getNameProps(obI.colStruc.names);


cellIDs = obI.nameProps.cellNum;
checkIDs = unique(cellIDs(cellIDs>0));
cellCount = 0;
for i = 1:length(checkIDs)
    obIDs = find((cellIDs == checkIDs(i)) & obI.nameProps.ofID);
    mainID = intersect(obIDs,find(obI.nameProps.cell));
    
    if length(mainID)>1
        obI.nameProps.names(mainID)
        mainID
        'too many main IDs'
        mainID = mainID(1);
    end
    
    if ~isempty(mainID)
        cellCount = cellCount+1;
        obI.cell.obIDs{cellCount} = obIDs;
        obI.cell.name(cellCount) = checkIDs(i);
        obI.cell.mainID(cellCount) = mainID;
        obI.cell.label{cellCount} = obI.nameProps.names{mainID};
        obI.cell.anchors(cellCount,:) =  obI.colStruc.anchors(mainID,:);
    end
    
end




save([MPN 'obI.mat'],'obI');