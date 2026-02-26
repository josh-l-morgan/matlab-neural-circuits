

    MPN = GetMyDir
colormap gray(256)
projectPieDir = [MPN 'morph\projectPie\'];
if ~(exist(projectPieDir,'dir')),mkdir(projectPieDir); end
    disp('loading')
    load([MPN 'obI.mat'])
    load([MPN 'dsObj.mat'])
    
    
    %%
    dim = 1;
    
    %%
    
    

%% Get Preferences
seedList = [108 201 170 110];
useList = obI2cellList_tracedTCR(obI);
conTo = makeConTo(obI,seedList);
seedPref = seedPreferences(seedList,useList);

cellList = useList.postList;
cellCon =  seedPref.cellPrefNoExclusion;
cellPref = seedPref.cellPrefNoExclusion(1,:)./sum(seedPref.cellPrefNoExclusion(1:2,:));
cellPref(isnan(cellPref)) = 0.5;
isPref = (seedPref.cellPrefNoExclusion(1,:)/max(seedPref.cellPrefNoExclusion(1,:)) + ...
seedPref.cellPrefNoExclusion(2,:)/max(seedPref.cellPrefNoExclusion(2,:)) );
    
 maxNum = 30;
    
    
    
    fsize = double(max(cat(1,dsObj.subs),[],1));

    
    %%
clf
scaleVec = 1./[4 4 4];
dims = [1 3];
imageRes = 6000;
newSize = ceil(fsize.* scaleVec);
newSize = newSize(dims);

minSub =  [Inf Inf];
maxSub = [0 0]
clear allSubs allCounts
for i = 1:length(cellList)
     targCell = cellList(i);
     
     
    cellSubs = getCellSubs(obI,dsObj,targCell);
    newSub = round(scaleSubs(cellSubs(:,dims),scaleVec(dims)));
    newInd = sub2ind(newSize(1:2),newSub(:,1),newSub(:,2));
    uind = unique(newInd);
    countInd = hist(newInd,uind);
    [y x] = ind2sub(newSize(1:2),uind);
    
    clear shiftSub
    shiftSub(:,1) = x - mean(x) + imageRes * isPref(i);
    shiftSub(:,2) = y - mean(y) + imageRes * cellPref(i);
    
    newMin = min(shiftSub,[],1);
    minSub = min(minSub, newMin);
    
    newMax = max(shiftSub,[],1);
    maxSub = max(newMax,maxSub);

     allSubs{i} = shiftSub;
     allCounts{i} = countInd;
     
end

     minSub = floor(minSub-1);
     fieldSize = ceil(maxSub + 1)-minSub;
     
     field = zeros(fieldSize);
 for i = 1:length(allSubs)
    
     shiftSubs = round(allSubs{i});
     shiftInds = sub2ind(fieldSize,shiftSubs(:,1)- minSub(1),shiftSubs(:,2)- minSub(2));
     field(shiftInds) = field(shiftInds) + [allCounts{i}]';
     
 end
     image(field)
     
  


