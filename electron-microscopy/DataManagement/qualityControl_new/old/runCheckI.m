clear all

TPN = GetMyDir

nams = GetPics(TPN)


for i = 1:length(nams)
    allI{i} = imread([TPN nams{i}]);
end

for i = 1:length(allI)
tiles(i) = checkIqual3s(allI{i});
end

%%
for i = 1:length(tiles)
    reportQ{i,1} = nams{i};
    reportQ{i,2} = tiles(i).quality(1);  
end

reportQ
