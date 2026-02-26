
TPN = GetMyDir

APN = findFolders(TPN);
allNams = {}; allTimes = []; allFolders = {};
for f = 1: length(APN)
    [aNams aTimes]= getTiles(APN{f});
    fID = ones(length(aNams),1) * f;
    allFolders = cat(1,allFolders,APN(fID));
    allNams = cat(1,allNams, aNams);
end

mkdir([TPN 'SubTiles'])

for i = 1:length(allNams)
   
useRegX = [4000 4300];
useRegY = [100 400];

I = 255-imread([allFolders{i} '\' allNams{i}],'PixelRegion',{useRegY,useRegX});
imwrite(I,[TPN 'SubTiles\subTile_' allNams{i}], 'Compression','none');


end






