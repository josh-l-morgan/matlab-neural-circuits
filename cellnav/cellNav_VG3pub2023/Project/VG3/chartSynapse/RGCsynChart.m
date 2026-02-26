


global tis

%% get list of vg3s
tis.cells.type.typeNames;
isAmc = (tis.cells.type.typeID == 8);
isVG3 = tis.cells.type.subTypeID == 1;
vg3sTarg = find(isAmc & isVG3);
vg3s = tis.cids(vg3sTarg)


%% count RGC syn
synAnchors = tis.obI.colStruc.anchors(tis.syn.obID,:);
for i = 1:length(vg3s)
   isPreTarg = tis.syn.pre == vg3s(i);
   synNum(i) = sum(isPreTarg);
   isRGC = tis.syn.postClass == 1;
   rgcSynNum(i) = sum(isPreTarg & isRGC);
   isOutUnk = tis.syn.postClass == 0;
   postUnkCids = tis.syn.post(isPreTarg & isUnk);
   unkNum(i) = length(postUnkCids);
   postRGCCids = tis.syn.post(isPreTarg & isRGC);
   uRGCs = unique(postRGCCids);
   rgcNum(i) = length(uRGCs);
   
   %% Get syn possition and ID
   L = unkNum(i);
   unkAnk = synAnchors(isPreTarg & isUnk,:);
   unkTypeCheckListCell{i} = cat(2,repmat(vg3s(i),[L 1]), postUnkCids, unkAnk);
   
    L = rgcSynNum(i);
   rgcAnk = synAnchors(isPreTarg & isRGC,:);
   rgcTypeCheckListCell{i} = cat(2,repmat(vg3s(i),[L 1]), postRGCCids, rgcAnk);
   
   
   
   
end



%% Format result

res(:,1) = vg3s;
res(:,2) = synNum;
res(:,3) = unkNum;
res(:,4) = rgcSynNum;
res(:,5) = rgcNum;

unkTypeCheckList = [];
rgcTypeCheckList = [];
for i = 1:length(vg3s)
   unkTypeCheckList = cat(1,unkTypeCheckList,unkTypeCheckListCell{i}) ;
   rgcTypeCheckList = cat(1,rgcTypeCheckList, rgcTypeCheckListCell{i});
   
end

