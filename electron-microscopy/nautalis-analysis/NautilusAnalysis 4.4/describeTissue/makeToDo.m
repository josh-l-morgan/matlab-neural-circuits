
load('MPN.mat')
%load('WPN.mat')

load([WPN 'tis.mat'])




isVG3idx = find((tis.cells.type.typeID==8) & (tis.cells.type.subTypeID == 1));
isVG3cids = tis.cells.cids(isVG3idx);

for i = 1:length(isVG3cids)
    
    synIn = find(tis.syn.post == isVG3cids(i));
    cellIn = setdiff(unique(tis.syn.pre(synIn)),0);
    
    synOut = find(tis.syn.pre == isVG3cids(i));
    cellOut = setdiff(unique(tis.syn.post(synOut)),0);
    
    
      
    
end



isBipIdx = find(tis.cells.type.typeID==7);
isBipCids = tis.cells.cids(isBipIdx)';
fvDir = [WPN 'fvLibrary\'];
clear faceNum
for i = 1:length(isBipCids)
   fileName = sprintf('%s\\%d.mat',fvDir,isBipCids(i));
   load(fileName) 
   faceNum(i,1) = size(fv.faces,1); 
end







