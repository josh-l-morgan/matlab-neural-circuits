function[cellSubs] = getCellSubs(obI,dsObj,cellId)

%%Compile all subs from each object id referenced by cellID

uCellId = unique(cellId(cellId>0));

cellSubs = [];
for i = 1:length(uCellId)
    targ = find(obI.cell.name==uCellId(i));
    obTarg = obI.cell.obIDs{targ};
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
        sub = double(dsObj(obTarg(o)).subs);
        if ~isempty(sub)
            cellSubs = cat(1,cellSubs,sub);
        end
        end
    end
end
