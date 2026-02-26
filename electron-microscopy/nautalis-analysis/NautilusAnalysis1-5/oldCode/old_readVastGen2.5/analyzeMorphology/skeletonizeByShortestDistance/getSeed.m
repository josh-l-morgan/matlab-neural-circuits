function[seed]   = getSeed(obI,cellName)


cellTarg = find(obI.cell.name == cellName);
obIDs = obI.cell.obIDs{cellTarg};

mainOb = obIDs(find(obI.nameProps.cell(obIDs),1));
seed = obI.colStruc.anchors(mainOb,:);