%clear all

utslDir = 'D:\MasterUTSL\UTSL_lgns1\';
cNum = 4;
rNum = 4;
%montageDir = 'Y:\joshm\LGNs1\rawMontages\'

%% Find all overviews
if ~exist('D:\MasterUTSL\UTSL_lgns1\tileList.mat');
    
    dirUTSL = dir(utslDir); dirUTSL = dirUTSL(3:end);
    wafNum = 0;
    
    for i = 1:length(dirUTSL)
        if dirUTSL(i).isdir
            nam = dirUTSL(i).name;
            disp(sprintf('Characterizing UTSL for %s',nam));
            if lower(nam(1)) == 'w'
                soDir = [utslDir nam '\SectionOverviewsDirectory\'];
                if exist(soDir,'dir') % is there an overview directory
                    dirSo = dir(soDir); dirSo = dirSo(3:end);
                    soNum = [];
                    wafNum = wafNum + 1;
                    wafSO(wafNum).name = nam;
                    
                    for s = 1:length(dirSo)  %find all section overviews
                        nam = dirSo(s).name;
                        und = regexp(nam,'_');
                        dotTif = regexp(nam,'.tif');
                        if ~isempty(und) & ~isempty(dotTif)
                            secNum = str2num(nam(und(1)+1:dotTif-1));
                            if ~isempty(secNum)
                                soNum(length(soNum)+1) = secNum;
                            end
                        end %if so file
                    end % run all files in so dir
                    
                    
                    wafSO(wafNum).soNum = soNum;
                    
                else  %if overview dir exists
                    
                end % if overview dir exists
            end % if wafer
        end % if dir
    end %end run dirUTSL
    save('..\matFiles\overviewList.mat','wafSO')

else
    load('..\matFiles\overviewList.mat')
end


%% Check for tiles
countS = 0;
clear utslTiles
for w = 1:length(wafSO)

    waf = wafSO(w).name;

    utslTiles.waf(w).name = waf;
    utslTiles.waf(w).overviews = wafSO(w).soNum;
    
    disp(sprintf('Checking for tiles in %s.',waf))
    for s = 1:length(wafSO(w).soNum)
        countS = countS+1
        sec = wafSO(w).soNum(s);
        
        secDir = sprintf('%s_Sec%03.0f_Montage',waf,sec); %[waf '_Sec' zeroBuf(sec) '_Montage'];
        
        
        utslTiles.waf(w).sec(s).secNum = sec;
        utslTiles.waf(w).sec(s).secDir = secDir;
        
        
        %%Generate tile Names
        tileNames = {};
        tifNames = tileNames;
        for r = 1:rNum
            for c = 1:cNum
                tileName = sprintf('Tile_r%d-c%d_%s_sec%03.0f', r, c, waf, sec);
                tileNames{length(tileNames)+1,1} = tileName;
                tifNames{length(tileNames),1} = [secDir '\' tileName '.tif'];
            end
        end
                   
        utslTiles.waf(w).sec(s).tileNames = tileNames;
        utslTiles.waf(w).sec(s).tifNames = tifNames;
        
    end % run all sections
    
end % run all wafers

save('..\matFiles\utslTiles.mat','utslTiles')


%% make single ID list
i = 0;
for w = 1:length(utslTiles.waf)
    for s = 1:length(utslTiles.waf(w).sec)
        for t = 1:length(utslTiles.waf(w).sec(s).tileNames)
            i = i+1;
            tileName = utslTiles.waf(w).sec(s).tileNames{t};
            utsltiles.tiles.tileName{i,1} = tileName;
            utslTiles.tiles.ID(i,1) = tile2ID(tileName); 
        end
    end
end

save('..\matFiles\utslTiles.mat','utslTiles')

