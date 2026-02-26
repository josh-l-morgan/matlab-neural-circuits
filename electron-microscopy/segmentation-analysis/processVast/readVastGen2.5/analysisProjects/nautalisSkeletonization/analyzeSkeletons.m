

MPN = GetMyDir
    aspect =  [1.2267 1.0667 1]

load([MPN 'obI.mat'])

skelDir = [MPN 'skel\']

morphDir = [MPN 'skel\morph\'];
if ~exist(morphDir,'dir'),mkdir(morphDir),end

skelMats = dir([skelDir 'mat\*.mat']);

%% pick Cells

useCells = obI.cell.name(obI.cell.isCell>0)
%tracedCells = []

%%
clear anaSkel
for i = 1:length(useCells)
    disp(sprintf('running cell %d of %d',i,length(useCells)))
    nam = sprintf('%d.mat',useCells(i));
    if exist([skelDir 'mat\' nam],'file')
        load([skelDir 'mat\' nam])
   
    cellID =useCells(i);
    seedSubA = cellStruct.surfVox.subs(cellStruct.surfVox.firstSeed,:)
    seedSubB = obI.cell.anchors(find(obI.cell.name == cellID),:)
    seedSub = seedSubB./[32 32 16];
    anaSkel(i).skelOri = plotSkeletons(cellStruct.arbor,seedSub);
    anaSkel(i).skelOri.cellID = cellID;
    pause(.01)
    anaSkel(i).skelOri
    anaSkel(i).cellID = cellID;
    
% 
%                 cellView = cellStruct.sideViews{1};
%                 skel = cellStruct.skel;
%                 
%                 
%                 save(skelFile,'cellStruct')
%                 %imwrite(cellView,viewFile)
%                 image(cellView)
%                 
%               
%         clf
fileName = sprintf('%smorphOri_%s.png',morphDir,nam(1:end-4))
saveas(gcf,fileName,'png')

 else
        ['cant find file ' nam]
    end


end
save([skelDir 'anaSkel.mat'],'anaSkel');


%% find traced cells
for i = 1:length(anaSkel)
    if isempty(anaSkel(i).skelOri)
        wasTraced(i) = 0;
    else
        
        foundTraced = find(tracedCells == anaSkel(i).skelOri.cellID,1);
        if isempty(foundTraced)
            wasTraced(i) = 0;
        else
            wasTraced(i) = 1;
        end
    end
    
end




%%
showTraced = find(wasTraced)
clear oriRats useLength
minorArbor = 1/3;
models = {[1 1 ; 1 1 ; 1 1] [minorArbor minorArbor ; 1 1; 1 1] [1 1; 1 minorArbor; 1 1] ...
    [1 1; minorArbor 1; 1 1]  [minorArbor minorArbor; minorArbor minorArbor; 1 1]}

for m = 1:length(models), models{m} = models{m}/sum(models{m}(:)),end
clear oriRats useLength tracedIDs
for i = 1:length(showTraced)
    
    
    skelOri = anaSkel(showTraced(i)).skelOri;
   oriRats(i,:) = skelOri.oriRat;
   useLength(i) = skelOri.lengthUsed;
   tracedIDs(i) = skelOri.cellID;
   
   
   dirCount = skelOri.dirCount;
   
   normCount = dirCount/sum(dirCount(:))
   
   for m = 1:length(models)
      difCount = models{m}-normCount
      modelDif(1,m) = sqrt(sum(sum(difCount.^2))) 
      
   end
    modelDifs(i,:)= modelDif;
   bestDif(i) = find(modelDif == min(modelDif));
   sortDif = sort(modelDif,'ascend');
   fitQual(i) = sortDif(2)-sortDif(1);
   subplot(2,1,1)
   totLength(i) = sum(dirCount(:));
   
    bar(normCount(:,1)*-1)
    hold on
    bar(normCount(:,2))
    hold off
    
    subplot(2,1,2)
      bar(modelDif)
      
   pause(.01)
   
   
end
oriRats(isnan(oriRats(:,1)),1) = 0;
oriRats(isnan(oriRats(:,2)),2) = 0;
showLength = useLength./max(useLength);

% 
% for i = 1:length(anaSkel.skelOri)
%     scatter(oriRats(i,1),oriRats(i,2),showLength(i));
%     hold on
% end
% hold off
clf
minLength = useLength>00;
scatter(oriRats(minLength,1),oriRats(minLength,2),showLength(minLength)*500+.001)
scatter(bestDif,fitQual)

reportFits = [uint16(tracedIDs') uint16(bestDif') uint16(fitQual' * 100) totLength']
