load('..\matFiles\utslTiles.mat')

%% make single ID list
i = 0;
for w = 1:length(utslTiles.waf)
    for s = 1:length(utslTiles.waf(w).sec)
        for t = 1:length(utslTiles.waf(w).sec(s).tileNames)
            i = i+1;
            tileName = utslTiles.waf(w).sec(s).tileNames{t};
            utslTiles.tiles.ID(i,1) = tile2ID(tileName); 
        end
    end
end