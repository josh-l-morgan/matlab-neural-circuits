function skelDat=getSkelDat(skelDir,cid)
    skelDat=struct;
    skelDat=load([skelDir 'sm_cid' num2str(cid) '.mat']);
end