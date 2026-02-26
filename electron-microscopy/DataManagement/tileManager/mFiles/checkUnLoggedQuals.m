
clear all

load('..\matFiles\logQuals.mat');
load('..\matFiles\tifMap.mat');


if exist('..\matFiles\noLogQuals.mat','file')
    load('..\matfiles\noLogQuals.mat')
else
    
    %% find un logged wafers
    
    
    noLogWafList = [];
    for w = 1:length(logQuals.noLogFound)
        nam = logQuals.noLogFound{w};
        noLogWafList(w) = str2num(nam(10:13));
    end
    
    noLogQuals.noLogWafList = noLogWafList;
    
    %% find unlogged tiles
    tileNum = 0;
    clear montageDir tifNames tileIDs sectionDirs subs
    for nL = 1:length(noLogQuals.noLogWafList);
        w = noLogQuals.noLogWafList(nL);
        for m = 1:length(tifMap.mon)
            IDs = tifMap.mon(m).tileIDs;
            monWafs = tifMap.mon(m).subs(:,1);
            targWaf = find(monWafs == w);
            for tW = 1:length(targWaf)
                t = targWaf(tW);
                tileNum = tileNum+1;
                noLogQuals.tiles(tileNum).montageDir = tifMap.mon(m).montageDir;
                noLogQuals.tiles(tileNum).tifNames = tifMap.mon(m).tifNames(t);
                noLogQuals.tiles(tileNum).tileIDs= tifMap.mon(m).tileIDs(t);
                noLogQuals.tiles(tileNum).sectionDirs = tifMap.mon(m).sectionDirs(t);
                noLogQuals.tiles(tileNum).subs = tifMap.mon(m).subs(t);
                noLogQuals.tiles(tileNum).checked = 0;
                
            end
        end
    end
    
    save('..\matFiles\.noLogQuals.mat','noLogQuals')
end

%%
%oldTileIDs = [noLogQuals.tiles.tileID];
for i = 1:length(noLogQuals.tiles);
    if ~noLogQuals.tiles(i).checked
        sprintf('Checking tile %d of %d.',i,length(noLogQuals.tiles))
        tileID = noLogQuals.tiles(i).tileIDs;
        tileName = ID2tile(tileID);
        tifNum = 0;
        for m = 1:length(tifMap.mon)
            tifPos = find(tifMap.mon(m).tileIDs == tileID);
            
            for t = 1:length(tifPos)
                
                targ = tifPos(t);
                tifPath = [tifMap.mon(m).montageDir '\' ...
                    tifMap.mon(m).sectionDirs{targ} '\' ...
                    tifMap.mon(m).tifNames{targ}];
            
                tifNum = tifNum + 1;
                noLogQuals.tiles(i).tif(tifNum).tifPath = tifPath;
                noLogQuals.tiles(i).tif(tifNum).tifName = tifMap.mon(m).tifNames{targ};
                noLogQuals.tiles(i).tif(tifNum).sectionDirs = tifMap.mon(m).sectionDirs{targ};
                noLogQuals.tiles(i).tif(tifNum).montageDir = tifMap.mon(m).montageDir;
                
                qual = checkFileQual(tifPath);
                noLogQuals.tiles(i).tif(tifNum).qual = qual;
                
            end  %check all tif positions
        end  %check all montage directories
        noLogQuals.tiles(i).checked = 1;
    end %if checked alread
    if mod(i,100)==0  %save periodically
        save('..\matFiles\.noLogQuals.mat','noLogQuals')

    end
end


save('..\matFiles\.noLogQuals.mat','noLogQuals')









