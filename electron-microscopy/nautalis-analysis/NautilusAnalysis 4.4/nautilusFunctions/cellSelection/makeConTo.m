function[conTo] = makeConTo(obI,targCells)

%%Targ cells is a vector of cell IDs mapped to obI

synapses = obI.nameProps.edges;
synapses = cat(1,[ 0 0 0], synapses);
allPostCells = unique(synapses(:,1));


for t = 1:length(targCells)
    conTo(t).targ = targCells(t);
end

tcrList = obI.nameProps.cellNum(obI.nameProps.tcr);
rgcList = obI.nameProps.cellNum(obI.nameProps.rgc);
linList = obI.nameProps.cellNum(obI.nameProps.lin);

for t = 1:length(conTo)
  
preWithTarg = unique(synapses(synapses(:,1)==conTo(t).targ,2));
preWithTarg = preWithTarg(preWithTarg>0);
rgcWithTarg = intersect(preWithTarg,rgcList);
linWithTarg = intersect(preWithTarg,linList);

synWithPre = [];
for i = 1:length(rgcWithTarg)
    synWithPre = cat(1,synWithPre, find((synapses(:,2)==rgcWithTarg(i)) & ...
        (synapses(:,1)>0)));
end
    
synPre = synapses(synWithPre,2);
synPost = synapses(synWithPre,1);
synObj = synapses(synWithPre,3);
postList = unique(synPost(synPost>0));

conTo(t).preList = preWithTarg;
conTo(t).rgcList = rgcWithTarg;
conTo(t).postList = postList;
conTo(t).tcrList = intersect(postList,tcrList);
conTo(t).postLIN = intersect(postList,linList);
conTo(t).preLIN = linWithTarg;

conTo(t).syn = [synPre synPost];
conTo(t).synObj = synObj;

end
