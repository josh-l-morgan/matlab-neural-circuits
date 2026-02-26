
TPN = GetMyDir;


    load([TPN 'q.mat'])
    
 APN = findFolders(TPN)
    allNams = {}; allTimes = []; allFolders = {};
for f = 1: length(APN)
    [aNams aTimes]= getTiles(APN{f});
    fID = ones(length(aNams),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
    allNams = cat(1,allNams, aNams);
    allTimes = cat(1,allTimes, aTimes);
end
   
    
    
checkit = ones(length(allNams),1);
for i = 1:length(allNams)
    for c = 1:length(q.checkedNams)
        if strcmp(allNams{i},q.checkedNams{c})
            q.checkedFolders{c} = allFolders{i};
            break
        end
    end
end


save([TPN 'q.mat'],'q');
