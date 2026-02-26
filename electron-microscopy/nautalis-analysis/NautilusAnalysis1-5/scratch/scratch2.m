

MPN = GetMyDir
load([MPN 'obI.mat'])

%% Find connectivity


seedList =[109 108 201 907 903 170];
conTo = makeConTo(obI,seedList);

[useList] = obI2cellList_seedInput(obI,seedList);
seedPref = seedPreferences(seedList, useList);

allEdges = obI.nameProps.edges(:,[2 1]);
cons = postTo(allEdges,9102)
cons = postTo(allEdges,9101)

%%

checkCell =232;
targ = find(seedPref.cellList == checkCell);
preCons = preTo(allEdges,checkCell)
conStats = [seedPref.sharedAx(:,targ) ... 
    seedPref.sharedSyn(:,targ)  ...
    seedPref.geoMeanNorm(:,targ)];

postTo(allEdges,checkCell);



%% find crossover
shared = {};
for x = 1:length(seedList)
    for y = 1:length(seedList)
        shared{y,x} = intersect(conTo(y).tcrList,conTo(x).tcrList);
    end
end

shared{2,3}
shared{2,4}
shared{2,5}
shared{3,4}
shared{3,5}

checkGiant = unique([shared{2,4} shared{2,5} shared{3,4} shared{3,5}])








