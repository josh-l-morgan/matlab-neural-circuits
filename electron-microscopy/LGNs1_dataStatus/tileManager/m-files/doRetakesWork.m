%%Compare qualities of multiply taken tiles to identify best.

SPN = 'D:\LGNs1\rawMontages_folder2\'; %source

dSPN = dir(SPN); dSPN = dSPN(3:end);

%% find wafers
c = 0;
clear wafNum sectionList
for i = 1:length(dSPN)
    nam = dSPN(i).name;
    secExp = regexp(nam,'_Sec');
    if secExp == 5
        c = c + 1;
        wafNum(c)= str2num(nam(2:4));
        sectionList{c} = nam;
    end
end

allWafs = sort(unique(wafNum),'descend')

%% run sections
for s = 1000:length(sectionList)
    
    dSec = dir([SPN sectionList{s} '\Tile*.tif']);
    tileNames = {dSec.name}';
    datenums = [dSec.datenum];
    [sortTimes tileOrder] = sort(datenums);
    tileList = tileNames(tileOrder);
    retook = 1:length(tileList)-1;
    for t = 1:length(tileList)-1
        tilePos = tileList{t}(1:10);
        retook(t) = strcmp(tilePos,tileList{t+1}(1:10));
    end
        disp(tileList)
        disp(retook)
        pause; 
    
    
end
   



%%


%%
missingTile = {};
for t = 1:length(gotWorse)
    tNam = gotWorse{t};
    slashes = regexp(tNam,'\');
    secNam = tNam(slashes(2)+1:slashes(3)-1);
    tileNam = tNam(slashes(3)+1:end-4);
    waf = str2num(secNam(2:4));
    if sum(find(allWafs == waf))
        
        
        if exist([SPN secNam],'dir')
            dSec = dir([SPN secNam])
            tilesToCheck = {};
            for f = 1:length(dSec)
                nam = dSec(f).name;
                if sum(regexp(nam,tileNam)) & sum(regexp(nam,'.tif'))
                    tilesToCheck{length(tilesToCheck)+1} = nam;
                end
            end
            
            if ~isempty(tilesToCheck)
                qualChecked = [];
                for tc = 1:length(tilesToCheck)
                    qualVal = checkFileQual([SPN secNam '\' tilesToCheck{tc}])
                    qualChecked(tc) = qualVal.quality;
                end
                bestQual{t} = tilesToCheck{find(qualChecked == max(qualChecked),1,'last')}
                
                
            else
                missingTile{length(missingTile)+1,1} = tileNam;
            end
            
            
        else
            
            
            
        end
    end
end

