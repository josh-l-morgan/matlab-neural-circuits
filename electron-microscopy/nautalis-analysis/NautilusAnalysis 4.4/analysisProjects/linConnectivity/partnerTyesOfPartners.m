


load('MPN.mat')
load([MPN 'obI.mat'])
mot = getMotifs(obI);
target = 125;


load([MPN 'nep\skelNep125.mat'])
load([MPN 'nep\useSynPos125.mat'])
pos = nep.nodePos;
skelEdges = nep.edges;

fromUnkPos = useSynPos.synPos125fromUNK;
edgeLength = getLengths(edges,pos);


fromUNKnum = size(useSynPos.fromUNK_dat,1);
fromRGCnum = size(useSynPos.fromRGC_dat,1);
fromLINnum = size(useSynPos.fromLIN_dat,1);
fromUNKnum/(fromRGCnum+fromLINnum)


%%

edges = obI.nameProps.edges(:,[1 2]);
t = mot.cel.types;

preTarg = unique(edges(find(edges(:,1) == target),2));
postTarg = unique(edges(find(edges(:,2) == target),1));

checkPre = intersect(preTarg,t.unks);
checkPre = checkPre((checkPre>9100)&(checkPre<10000));

clear postOfPre postOfPreType
for i = 1:length(checkPre)
    
    parts = setdiff(unique(edges(find(edges(:,2) == checkPre(i)),1)),[target 0]);
    
    partsType = [];
    for p = 1 : length(parts)
        partsType(p,:) = [sum(t.rgcs == parts(p)) sum(t.tcrs == parts(p)) ...
            sum(t.lins == parts(p)) sum(t.unks == parts(p))];
    end
    postOfPre{i} = parts(:);
    postOfPreType{i} = partsType;
    
end

partnerTypes = cat(1,postOfPreType{:});
partners = cat(1,postOfPre{:});
report = [partners partnerTypes]
sum(partnerTypes)

scatter3(pos(:,1),pos(:,2),pos(:,3),'.')
hold on
scatter3(fromUnkPos(:,1),fromUnkPos(:,2),fromUnkPos(:,3),'o','filled','r')

hold off

