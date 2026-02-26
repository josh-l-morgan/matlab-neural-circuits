

cladeRes = result;

%% Get groups
groupNum = 2;
groupNodes = getSplit(cladeRes,groupNum)




%% Get recorded properties
[checkIDs checkProp]  = getList_axSeedCon();





[checkIDs checkProp] = getList_giantBoutonsCounts(MPN);
useIDs = checkIDs<1000;
checkIDs = checkIDs(useIDs);
checkProp = checkProp(useIDs);
checkProp(checkProp>0) = 1;% checkProp(checkProp>0) + max(checkProp);



[checkIDs checkProp] = getList_pcLatRats;
checkProp = checkProp(:,1);



[groupPropIDs groupProps] = getGroupProps(groupNodes,checkIDs,checkProp);

%%
clf
allProps = cat(1,groupProps{:});


groupMat = [];
for i = 1:length(groupProps)
   props = groupProps{i};
   groupMat = cat(1,groupMat,[repmat(i,[size(props,1) 1]) props]);    
end


for i = 1:length(groupProps)
    props = groupProps{i};
    gs(i).N = size(props,1);
    gs(i).dim = size(props,2);
    gs(i).mean = mean(props,1);
    gs(i).std = std(props,1);
    gs(i).stderror = std(props,1)/sqrt(size(props,1));
    gs(i).variance = var(props,0,1);
end

gComp.gs = gs
gComp.ranksum = ranksum(groupProps{1},groupProps{2})
%gComp.KW = kruskalwallis(groupMat)

binProp = [min(allProps):range(allProps)/10:max(allProps)];



hist1 = hist(groupProps{1},binProp);
hist2 = hist(groupProps{2},binProp);
subplot(1,1,1)
bar(binProp,[hist1; hist2]')
ranksum(groupProps{1},groupProps{2})


