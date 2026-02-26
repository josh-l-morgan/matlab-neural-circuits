function [x1,y1,z1,x2,y2,z2]=getBoundingBox(cid)
global tis
if isempty(tis)
    'This requires a tis.mat to be loaded'
else
    allObjCids=tis.obI.nameProps.cellNum;
    cidObjIDs=find(allObjCids==cid);
    cidBoxes=tis.obI.colStruc.boundBox(cidObjIDs,:);
    
end