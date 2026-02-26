
clear all

load('..\matFiles\logQuals.mat');
load('..\matFiles\tifMap.mat');



%% Collect suspect tiles

if exist('..\matFiles\chkQ.mat')
    load('..\matFiles\chkQ.mat')
else

%if ~isfield(chkQ,'do')
worseThresh = .2;
countTiles = 0;
clear tileNames tileIDs tileWorse
for w = 1: length(logQuals.waf)
    if logQuals.waf(w).foundBook
        for s = 1:length(logQuals.waf(w).sec)
            for t = 1:length(logQuals.waf(w).sec(s).tile)
                
                countTiles = countTiles + 1;
                tileName = logQuals.waf(w).sec(s).tile(t).tileName;
                tileNames{countTiles,1} = tileName;
                tileIDs(countTiles,1) = tile2ID(tileName);
                tileWorse(countTiles,1) = (logQuals.waf(w).sec(s).lastQuals(t) + ...
                    logQuals.waf(w).sec(s).lastQuals(t) * worseThresh) < ...
                    logQuals.waf(w).sec(s).bestQuals(t);
            end
        end
    end
    %end
    
    chkQ.logRes.tileNames = tileNames;
    chkQ.logRes.tileIDs = tileIDs;
    chkQ.logRes.tileWorse = tileWorse;
    
    chkQ.do = tileIDs(tileWorse>0);
    chkQ.done = zeros(length(chkQ.do),1);
    
    save('..\matFiles\chkQ.mat','chkQ')
end
end

%%

%oldTileIDs = [chkQ.tiles.tileID];
for i = 1:length(chkQ.do);
    sprintf('Checking tile %d of %d.',i,length(chkQ.do))
    tileID = chkQ.do(i);
    tileName = ID2tile(tileID)
    chkQ.tiles(i).tileID = tileID;
    chkQ.tiles(i).tileName = tileName;
    tifNum = 0;
    for m = 1:length(tifMap.mon)
        tifPos = find(tifMap.mon(m).tileIDs == tileID);
        
        for t = 1:length(tifPos)
            
            
            targ = tifPos(t);
            tifPath = [tifMap.mon(m).montageDir '\' ...
                tifMap.mon(m).sectionDirs{targ} '\' ...
                tifMap.mon(m).tifNames{targ}];
            
            % Check to see if the tif has been checked before
%             oldPaths = chkQ.tiles(i).tif(tifNum).tifPath;
%             alreadyChecked = 0;
%             for p = 1:length(oldPaths)
%                 alreadyChecked = strcmp(oldPaths{p},tifPath);
%                 if alreadyChecked
%                     break
%                 end
%             end
            
            
            tifNum = tifNum + 1;
            chkQ.tiles(i).tif(tifNum).tifPath = tifPath;
            chkQ.tiles(i).tif(tifNum).tifName = tifMap.mon(m).tifNames{targ};
            chkQ.tiles(i).tif(tifNum).sectionDirs = tifMap.mon(m).sectionDirs{targ};
            chkQ.tiles(i).tif(tifNum).montageDir = tifMap.mon(m).montageDir;
            
            qual = checkFileQual(tifPath);
            chkQ.tiles(i).tif(tifNum).qual = qual;
            
        end
        
    end
    
    chkQ.done(i) = 1;
end



save('..\matFiles\chkQ.mat','chkQ')

