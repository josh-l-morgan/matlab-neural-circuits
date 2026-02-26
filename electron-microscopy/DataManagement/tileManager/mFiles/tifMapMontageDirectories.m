%% Enter data locations

load('..\matFiles\montageDirs.mat')
if exist('..\matFiles\tifMap.mat','file')
    load('..\matFiles\tifMap.mat')
end
for m = 2:length(montageDirs)
           montageDir = montageDirs{m};
            mDir = dir(montageDir); 
            tifCount = 0;
            tifNames = cell(200000,1);
            tifTimes = zeros(200000,1);
            tileNames = cell(200000,1);
            tileIDs = zeros(200000,1);
            sectionDirs = cell(200000,1);
            
            tifMap.mon(m).montageDir = montageDir;
            
            for s = 1:length(mDir)  % run all sections
                nam = mDir(s).name;
                sprintf('mapping section %d of %d',s,length(mDir))
                if sum(regexp(nam,'_Montage'))
                    sDir = dir([montageDir '\' nam]);
                    names = cat(1,{sDir.name});
                    dates = datenum(cat(1,sDir.date));
                    tiles = ~cellfun(@isempty,regexp(names,'Tile_'));
                    tifs = ~cellfun(@isempty,regexp(names,'.tif'));
                    checkTiles = find(tiles & tifs);
                    
                    
                    for t = 1:length(checkTiles)
                       tifCount = tifCount+1; 
                       tifName = names{checkTiles(t)}; 
                       secPos = regexp(tifName,'_sec');
                       tileName = tifName(1:secPos+6);
                       
                       tifNames{tifCount} = tifName;
                       tifTimes(tifCount) = dates(checkTiles(t));
                       tileNames{tifCount} = tileName;
                       tileIDs(tifCount) = tile2ID(tileName);
                       sectionDirs{tifCount} = nam;
                    end
                    
                end
            end
            
            tifMap.mon(m).tifNames = tifNames(1:tifCount);
            tifMap.mon(m).tifTimes = tifTimes(1:tifCount);
            tifMap.mon(m).tileNames = tileNames(1:tifCount);
            tifMap.mon(m).tileIDs = tileIDs(1:tifCount);
            tifMap.mon(m).sectionDirs = sectionDirs(1:tifCount);
            
end

save('..\matFiles\tifMap.mat','tifMap')