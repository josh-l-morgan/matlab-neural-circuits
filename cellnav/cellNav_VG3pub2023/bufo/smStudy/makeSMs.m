

global tis

isAMC = find(tis.cells.type.typeID == 8); 
isVGC = intersect(isAMC,find(tis.cells.type.subTypeID == 1));
vgcCids = tis.cids(isVGC)

makeVolMPNcnv
smCells = vgcCids;
try close(smFig),end
smFig = figure
for i = 2:length(smCells)
    
    makeSM(smCells(i));
end


makeSM(1008)