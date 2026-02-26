function[boutList] = useListBoutSize(useList);


%%
load('MPN.mat');
load([MPN 'butSize4.mat']);
load([MPN 'obI.mat']);



axList = butSize.axList;
edges = butSize.edges;
butSize.synDat.edges;



%%
seedList = useList.seedList;
seedNum = length(seedList);
%useList = obI2cellList_seedInput_RGC_TCR(obI,seedList);
conTo = makeConTo(obI,seedList);

postList = useList.postList;
edges = butSize.edges;
con = edge2con(edges);
axList = butSize.axList;

%{
conPref =seedPreferences(seedList,useList);
isMixed = sum(conPref.sharedSyn([1 2],:)>0,1)==2;
mixList = conPref.cellList(isMixed);
%mixList = postList(randperm(length(postList),26));
unMixList = setdiff(postList,mixList);
%}

% %% Get size matrix
% clear compSizes
% for a = 1:length(axList)
%    compSizes(a,:) = cat(2,a, size(butSize.butVols{a},2),...
%    size(butSize.axProfile{a},1),...
%    size(butSize.synDat(a).edges,1));
% end

clear preSeed onSeedVols offSeedVols onMixVols
boutCon = useList.con * 0;
for a = 1:length(axList) % check each axon
        
    axID = axList(a);
    butVols = butSize.butVols{a};
    
    bSize = (butVols*3/4/pi).^(1/3)*butSize.voxLength * 2;
    bSize = pi * bSize.^2;
    %bSize = butVols * butSize.voxVol;
    bSize = butSize.butDiam{a};
    
    postBut = butSize.synDat(a).edges(:,1);
    
    preListTarg = find(useList.preList == axID);
    
    if ~isempty(preListTarg)
        
        for p = 1:length(postBut)
            postListTarg = find(useList.postList == postBut(p));
            if ~isempty(postListTarg)
                
                boutCon(preListTarg,postListTarg) = boutCon(preListTarg,postListTarg) + bSize(p);
                
            end
        end
    end
    
end

image(boutCon*20)

boutList.seedList = useList.seedList;
boutList.preList = useList.preList;
boutList.postList = useList.postList;
boutList.con = boutCon;



