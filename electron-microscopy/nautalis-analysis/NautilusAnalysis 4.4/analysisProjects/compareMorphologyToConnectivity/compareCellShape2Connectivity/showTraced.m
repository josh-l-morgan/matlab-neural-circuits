


    MPN = GetMyDir

projectPieDir = [MPN 'morph\projectPie2\'];
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
isPref = sum(seedPref.cellPrefNoExclusion(1:2,:));
    
 maxNum = 30;
    
    
    
    
    
    %%
clf
for i = 1:length(cellList)
     targCell = cellList(i);
    showCellNames = {targCell}  
    col = [1 1 1];   
 
 

fsize = double(max(cat(1,dsObj.subs),[],1));
minVal = double(min(cat(1,dsObj.subs),[],1));
viewProps.viewWindow = [minVal; fsize];

viewProps.dim = dim;
viewProps.maxScaleFactor = .1;
viewProps.sumScaleFactor = 3;
viewProps.obI = obI;
viewProps.dsObj = dsObj;
viewProps.col = col;
viewProps.fsize = fsize;
viewProps.cellId = showCellNames;









%%Display Cells

I= showCellsAndMore(viewProps);
[y x] = find(sum(I,3)>0);
I = I(min(y):max(y),min(x):max(x),:);
meany = mean(y)-min(y);
meanx = mean(x) - min(x);


image(uint8(256-I*1))
clear p
vals = cellCon(1:2,i);
if sum(vals)>= maxNum
    vals = [vals * maxNum/sum(vals); 0];
else
    vals = [vals; maxNum - sum(vals)];
end
p.vals = vals
%p.vals = [10 3 1];
p.color = [1 0 1; 0 1 0; 0 0 0]/3;
p.center = [size(I,1)/2 size(I,2)*2];
p.center = [meany meanx];
p.radius = 300;

hold on
plot([])
myPie(p)
hold off
pause(1)
% 
% imshow(I)
% hold on
% h = imshow(I);
% alpha = I>0;
% set(h, 'AlphaData', alpha);
% 
% 
%  imshow(E, 'InitialMag', 'fit')



imageName = sprintf('%ssumPie_%03.0f.eps',projectPieDir,targCell)
saveas(gcf,imageName,'eps')

end
%%










