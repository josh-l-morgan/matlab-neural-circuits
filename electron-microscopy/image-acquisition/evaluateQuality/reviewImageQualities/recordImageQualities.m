

%% Get source directory
SPN = GetMyDir;
dSPN = dir(SPN); dSPN = dSPN(3:end);

%% Find sections in source directory
secList = {};
for i = 1:length(dSPN);
    nam = dSPN(i).name;
    if sum(regexp(nam,'_Montage')) && dSPN(i).isdir
        secList{length(secList)+1} = nam;
    end
end


%%
allFlips = {};
for i = 1:length(secList)
    secDir = [SPN secList{i}]
    dSec = dir(secDir); dSec = dSec(3:end);
    if ~exist([secDir '\bestTile.mat'],'var')
        tileList = {};baseList = {};
        for t = 1:length(dSec)
            nam = dSec(t).name;
            if sum(regexp(nam,'Tile_r')) && sum(regexp(nam,'.tif'))
                tileList{length(tileList)+1,1} = nam;
                baseList{length(tileList),1} = nam(1:22);
            end
        end
        
        [uniqueBase ia ic] = unique(baseList); %find unique bases and their possitions
        histBase = hist(ic,[1:1:max(ic)]);
        bestTile  = cell(4,4);
        for b = 1:length(uniqueBase)
            nam = uniqueBase{b};
            underPos = regexp(nam,'_');
            rPos = regexp(nam,'_r');
            cPos = regexp(nam,'-c');
            r = str2num(nam(rPos+2:cPos-1));
            c = str2num(nam(cPos+2:underPos(2)-1));
            
            if (r*c)>0 %if this is a tile position
                checkTiles = find(ic==b);
                if length(checkTiles)==1
                    bestTile{r,c} = tileList{checkTiles};
                else
                    tileQuals = zeros(1,length(checkTiles));
                    for cT = 1:length(checkTiles)
                        %[secDir '\' tileList{checkTiles(cT)}]
                        qual = checkFileQual([secDir '\' tileList{checkTiles(cT)}]);
                        tileQuals(cT) = qual.quality;
                    end
                    qTarg = find(tileQuals == max(tileQuals),1,'last');
                    bestTile{r,c} = tileList{checkTiles(qTarg)};
                    if length(tileList{checkTiles(qTarg)})>25 %record bad retakes
                        allFlips{length(allFlips)+1} = tileList{checkTiles(qTarg)};
                    end
                end
                
                
            end %if tile position
          
            
        end 
          bestTile
            save([secDir '\bestTile.mat'],'bestTile')
        
    end%if bestTile does not exist
    %
    %    for t = 1:length(tileList)
    %       qual = checkFileQual([secDir '\' tileList{t}])
    %    end
    %
    
    
    
end

