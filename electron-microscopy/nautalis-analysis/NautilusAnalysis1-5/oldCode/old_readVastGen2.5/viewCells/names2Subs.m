function[cellSubs] = names2Subs(obI,dsObj,cellId);


%%
uCellId = cellId;

for i = 1:length(uCellId)
    
    term = uCellId{i};
    if strcmp(class(term),'char')
        for n = 1:length(obI.nameProps.names)
            isTerm(n) = sum(regexp(lower(obI.nameProps.names{n}),lower(term)));
        end
        obTarg = find(isTerm);
    else
        targ = find(obI.cell.name==term);
        obTarg = obI.cell.obIDs{targ};
    end
    stackSubs = [];
    for o = 1:length(obTarg)
        if obTarg(o)<=length(dsObj)
            sub = double(dsObj(obTarg(o)).subs);
            
            if ~isempty(sub)
                
                stackSubs = cat(1,stackSubs, sub);
            end
            
        end
        
    end
    
    cellSubs{i} = stackSubs;
end
