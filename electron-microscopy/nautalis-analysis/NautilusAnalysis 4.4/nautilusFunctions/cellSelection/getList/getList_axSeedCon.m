function[checkIDs checkCol]  = getList_axSeedCon();

MPN ='D:\LGNs1\Analysis\segmentations\vastMip3\collectVastSubs_mat\';

load([MPN 'obI.mat'])
seedList = [ 108 201 109 907 903];
useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
synEdges = obI.nameProps.edges(:,[2 1]);

checkIDs = useList.preList;

seedCol = [1 0 0; 0 1 0;  0 0 1; 0 .0 1; 0 0 1];
checkCol = zeros(length(checkIDs),3);
for s = 1:length(seedList)
    useCol = seedCol(s,:);
    for n = 1:length(checkIDs)
        checkCol(n,:) = checkCol(n,:) + useCol *  ...
            double(sum((synEdges(:,1)==checkIDs(n)) & (synEdges(:,2)) == seedList(s))>0);
    end
end
