function[obI] = fixColStruc(obI);

global glob


%%Embarassing hak
if strcmp(glob.NA.export.exportName,'Josh_B') & ...
    strcmp(glob.vol.activeName,'v1') 
    hackCellNavGlob
end

if isfield(glob.NA.export,'fix')

    %%Shift z position of all subs
    if isfield(glob.NA.export.fix,'changeZ')
        obI.colStruc.anchors(:,3) = obI.colStruc.anchors(:,3) + glob.NA.export.fix.changeZ;
        obI.colStruc.anchors = max(obI.colStruc.anchors,-1);
    end

end


