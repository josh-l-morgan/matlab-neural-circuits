function[cellNames cellProps] = getList_ldm(obI)


allIDs = [obI.nameProps.allIDs{obI.nameProps.ldm}];

cellNames = unique(allIDs)';
cellProps = hist(allIDs,cellNames)';
