function newFVs=flipFV(fvs)
if isfield(fvs,'vertices')
    fvs.vertices=fvs.vertices(:,[3 1 2]);
    newFVs=fvs;
else
    for i=1:size(fvs,1)
        curFV=fvs{i};
        curFV.vertices=curFV.vertices(:,[3 1 2]);
        newFVs{i}=curFV;
    end
end
end