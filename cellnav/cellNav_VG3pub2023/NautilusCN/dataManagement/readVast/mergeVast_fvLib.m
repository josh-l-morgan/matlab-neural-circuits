function[obI] = mergeVast_fvLib(mergeDirs,fvDir,vols)


         

%% parse names
fuse.exportDir = mergeDirs;



%% merge subs
lastInd = 0;
%allSubs = {};
clear allOb
allOb.fuse = fuse;
obSource = [];
for i = 1:length(fuse.exportDir)
    %load(mergeDirs{i} '\vastSubs.mat']);
    load([mergeDirs{i} '\obI.mat']);

    numOb = length(obI.colStruc.names);
    
    obSource(length(obSource)+1:length(obSource)+numOb) = i;
    if i == 1
        allOb.colStruc.info = obI.colStruc.info;
        allOb.colStruc.names = obI.colStruc.names;
        allOb.colStruc.ids = obI.colStruc.ids;
        allOb.colStruc.anchors = obI.colStruc.anchors;
        allOb.colStruc.boundBox =  obI.colStruc.boundBox;
        allOb.colStruc.col1 = obI.colStruc.col1;
        allOb.colStruc.col2 = obI.colStruc.col2;
        allOb.colStruc.parent = obI.colStruc.parent;
        allOb.colStruc.children = obI.colStruc.children;
        %allSubs = vastSubs(1:numOb);
        allOb.em = obI.em;
    else
        
        allOb.colStruc.names = cat(2,allOb.colStruc.names,obI.colStruc.names);
        ids = obI.colStruc.ids + lastInd;
        allOb.colStruc.ids = cat(2,allOb.colStruc.ids,ids);
        allOb.colStruc.anchors = cat(1,allOb.colStruc.anchors,obI.colStruc.anchors);
        allOb.colStruc.boundBox = cat(1,   allOb.colStruc.boundBox, obI.colStruc.boundBox);
        allOb.colStruc.col1 = cat(1,allOb.colStruc.col1,obI.colStruc.col1);
        allOb.colStruc.col2 = cat(1,allOb.colStruc.col2,obI.colStruc.col2);
        allOb.colStruc.parent = cat(2,allOb.colStruc.parent,obI.colStruc.parent);
        allOb.colStruc.children = cat(2,allOb.colStruc.children,obI.colStruc.children);
        
        %allSubs = cat(2,allSubs,vastSubs(1:numOb));
    end
    lastInd = max(allOb.colStruc.ids);
    allOb.fuse.mip(i) = obI.em.mipLevel;
    allOb.mergeFv(i).em = obI.em;
    allOb.mergeFv(i).fuse = obI.fuse;
    allOb.mergeFv(i).source = mergeDirs{i};
    allOb.mergeFv(i).volName = vols{i};
end

allOb.fuse.obSource = obSource;
obI = allOb;
%vastSubs = allSubs;

%%

%% choose dat
aliases = {};
for i = 1:length(mergeDirs)
    load([mergeDirs{i} '\tis.mat']);
    try
        aliases = {aliases{:} tis.dat.alias{:}};
    end
end


obI.nameProps = getNameProps2019(obI.colStruc.names,aliases);
obI = getObiCellProps(obI);


%save([MPN 'vastSubs.mat'],'vastSubs','-v7.3')
save([fvDir 'obI.mat'],'obI')

%% Down sample vastSubs
% 
% 
% mipLevs = obI.fuse.mip(obI.fuse.obSource);
% 
% res = obI.em.res; 
% vRes = [2.^mipLevs' 2.^mipLevs' mipLevs'*0+1].* res;
% dsRes = [100 100 100];
% dsDim = repmat(dsRes,[size(vRes,1) 1])./vRes;
% 
% dsObj = downSampObj(MPN, dsDim);
% % end
% 
% 
% obI.em.res = res; %or [4.6 4 30]?
% obI.em.vRes =vRes;
% obI.em.dsRes =dsRes/1000;
% save([MPN 'obI.mat'],'obI')




