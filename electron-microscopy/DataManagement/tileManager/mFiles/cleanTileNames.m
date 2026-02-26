
load('..\matFiles\tifMap.mat');


for m = 1:length(tifMap.mon)
    clear tileNames tileIDs subs
    tifNames = tifMap.mon(m).tifNames;
    for t = 1:length(tifNames)
        tileNames{t,1} = tif2tile(tifNames{t});
        [ID w s r c] = tile2ID(tileNames{t});
        tileIDs(t,1) = ID;
        subs(t,:) = [w s r c];
    end
    tifMap.mon(m).tileNames = tileNames;
    tifMap.mon(m).tileIDs = tileIDs;
    tifMap.mon(m).IDchar = num2str(tileIDs);
    tifMap.mon(m).subs = subs;
end