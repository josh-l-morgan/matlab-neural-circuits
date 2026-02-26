function[soma] = obIobSoma(obI,dsObj,cid)




targ = obI.cell.name == cid;
soma.obIDs = intersect(obI.cell.obIDs{targ}, find(obI.nameProps.tag.soma));

if ~isempty(soma.obIDs)
    
    soma.subs = double(cat(1,dsObj(soma.obIDs).subs));
    if isempty(soma.subs)
        soma.center = [];
    else
        soma.center = mean(soma.subs,1);
    end
    
else
    
    cellTarg = find(obI.cell.name == cid);
    mainID = obI.cell.mainID(cellTarg);
    anch = double(obI.cell.anchors(cellTarg,:));
    soma.center = anch;
    soma.subs = anch;
end

if sum(soma.center>0) ~= 3
    
    ids =  obI.cell.obIDs{targ};
    
    subs = double(cat(1,dsObj(ids).subs));
    
    
    if ~isempty(subs)
        modeS = mode(subs,1);
        
        dist = sqrt((subs(:,1)-modeS(1)).^2 + ...
            (subs(:,2)-modeS(2)).^2 + ...
            (subs(:,3)-modeS(3)).^2);
        
        closest = find(dist == min(dist),1);
        soma.center = subs(closest,:);
        soma.subs = subs(closest,:);
    else
        soma.center = [1 1 1];
        soma.subs = [1 1 1];
    end
    
end









