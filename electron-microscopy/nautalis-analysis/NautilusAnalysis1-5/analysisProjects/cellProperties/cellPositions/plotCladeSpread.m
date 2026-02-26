         




clear all
load('MPN.mat')
load([MPN 'obI.mat'])

[cellIDs cellPos] = parseCellPositions(obI)

load('.\data\clade\cladePick_six2.mat')

targCells = [108 201 903 907];
conTo = makeConTo(obI,targCells);
useCells = [conTo([1:4]).tcrList];


colorTable = [ 1 1 0; 1 0 1; 1 0 0; 0 .8 .8; .2 .2 1 ; 0 1 0] 
colorTable = [ .4 .3 1; 1 0 1; 1 0 0; 0 .8 .8; .2 .2 1 ; 0 1 0] 
colorTable = [ .2 .2 1; 1 1 0; 1 0 0; .5 .5 .5; 1 0 1 ; 0 1 0]

col = [];
allMembers = [];
useClade = [ 1 2 3 4 5 6];
for i = 1:length(useClade)
   members = cladePick.members{useClade(i)};
   members = intersect(useCells,members);
   allMembers = [allMembers members];
   col = cat(1,col,repmat(colorTable(useClade(i),:),[length(members),1]));
   
   tempPos = [];
   for m = 1:length(members)
      targ = find(cellIDs == members(m));
      tempPos(m,:) = cellPos(targ,:);
   end
   memPos{i} = tempPos;
   
end

%%
clf
hold on
histBin = [0:20:400];
clear stdev
for i = 1:length(memPos)
    pos = memPos{i};
    meanPos = mean(pos,1);
    dists = sqrt((meanPos(:,1)-pos(:,1)).^2 + ...
        (meanPos(:,2)-pos(:,2)).^2 + ...
        (meanPos(:,3)-pos(:,3)).^2);
    histDist = histc(dists,histBin);
    plot(histDist,'color',colorTable(i,:))
    stdev(i,:) = std(pos,1)
    memMean(i) = mean(dists);
end
hold off


%%
clf
hold on
histBin = [0:20:400];
clear stdev
for i = 1:length(memPos)
    pos = memPos{i};
    meanPos = mean(pos,1);
    dists = sqrt((meanPos(:,1)-pos(:,1)).^2 + ...
        (meanPos(:,2)-pos(:,2)).^2);% + ...
        %(meanPos(:,3)-pos(:,3)).^2);
    histDist = histc(dists,histBin);
    plot(histDist,'color',colorTable(i,:))
    stdev(i,:) = std(pos,1)
end
hold off


%%
clf
hold on
histBin = [0:20:400];
for i = 1:length(memPos)
    pos = memPos{i};
    meanPos = mean(pos,1);
    dists = sqrt((meanPos(:,2)-pos(:,2)).^2);% + ...
        %(meanPos(:,2)-pos(:,2)).^2);% + ...
        %(meanPos(:,3)-pos(:,3)).^2);
    histDist = histc(dists,histBin);
    plot(histDist,'color',colorTable(i,:))
    
end
hold off

