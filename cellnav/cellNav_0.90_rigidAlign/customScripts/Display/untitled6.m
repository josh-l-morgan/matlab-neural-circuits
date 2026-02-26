%%Apply segmentation alignment



obI.nameProps = getNameProps2019(obI.colStruc.names);

obI = getObiCellProps(obI);

save([MPN 'vastSubs.mat'],'vastSubs','-v7.3')
save([MPN 'obI.mat'],'obI')

%% Down sample vastSubs


mipLevs = obI.fuse.mip(obI.fuse.obSource);

res = obI.em.res; 
vRes = [2.^mipLevs' 2.^mipLevs' mipLevs'*0+1].* res;

try dsRes = glob.NA.export.dsRes*1000;
catch err
    dsRes = 0.1;
end
dsRes = max(dsRes,max(vRes(:)));
dsDim = repmat(dsRes,[size(vRes,1) 1])./vRes;

dsObj = downSampObj(MPN, dsDim);
% end


obI.em.res = res; %or [4.6 4 30]?
obI.em.vRes =vRes;
obI.em.dsRes =dsRes/1000;
save([MPN 'obI.mat'],'obI')


%% Clean synapses
[badSegList,overlapDist,newObi,sourceDirs]=cleanObi(obI,dsObj);

badSegList = [0; badSegList];
keepSynProp = ones(length(obI.nameProps.synProp),1);
for i = 1:length(obI.nameProps.synProp)
    try segID = obI.nameProps.synProp{i}.segID;
    catch err 
        segID = 0;
    end

    if sum(badSegList==segID)
        keepSynProp(i)= 0;
    end

end
obI.nameProps.synProp = obI.nameProps.synProp(keepSynProp>0);


edges = obI.nameProps.edges;
keepEdges = ones(size(edges,1),1);

for i = 1:size(edges,1)
    if sum(badSegList==edges(i,3))
        keepEdges(i) = 0;
    end
end
edges = edges(keepEdges>0,:);
obI.nameProps.edges = edges;


save([MPN 'obI.mat'],'obI')


