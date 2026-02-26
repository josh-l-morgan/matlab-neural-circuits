global tis






tis.cells.type.typeID == 1;


cidsPostToVGC = unique(COI.postToVGC);
COI.rgcGroupCids
rgcNames = {COI.rgcGroupNames{1}}
COI.typesPostToVGC


numCids = [];
subCids = [];
for i = 1:length(COI.rgcGroupCids)
    numCids(i) = length(COI.rgcGroupCids{i});
    subCids = cat(1,subCids,COI.rgcGroupCids{i});
end
totCids = sum(numCids)

subCidPostVGC = intersect(subCids,cidsPostToVGC);