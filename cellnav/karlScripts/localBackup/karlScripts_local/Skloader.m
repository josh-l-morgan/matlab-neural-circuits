%Just the facts, maam
vgcSkelCidList=[2 3 4 5 13 14];
skelDir='Y:\karlsRetina\CellNavLibrary_IxQ\Volumes\AprilMerge\Analysis\SMs\';
allSMs={};
for i=1:length(vgcSkelCidList)
    fileName=['sm_cid' num2str(vgcSkelCidList(i)) '.mat'];
    curSM=load([skelDir fileName]);
    allSMs{i}=curSM(:);
end

exportStruct={};
for i=1:length(vgcSkelCidList)
    curSkel=allSMs{i};