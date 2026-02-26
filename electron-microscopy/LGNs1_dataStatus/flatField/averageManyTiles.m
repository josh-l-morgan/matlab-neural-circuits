
SPN = 'D:\Other\flatFieldWaf1\';

%% Find tiles
dSPN = dir(SPN); dSPN = dSPN(3:end);    

tileNames = {};
foldNames = {};
for i = 1:length(dSPN)
    if dSPN(i).isdir
        foldNam = dSPN(i).name
        if regexp(foldNam,'Montage')
            dFold  = dir([SPN foldNam]); dFold = dFold(3:end);
            for f = 1:length(dFold)
                nam = dFold(f).name;
                if regexp(nam, 'Tile') & regexp(nam,'.tif')
                   
                    tileNames{length(tileNames)+1,1} = nam;
                     foldNames{length(tileNames),1} = [foldNam '\'];
                end
            end
        end
    end
end


%% Read tiles
numTiles = length(tileNames);

clear allI
for i = 1:numTiles
    disp(sprintf('averaging %d of %d',i,numTiles))
    %I = double(imread([SPN foldNames{i} tileNames{i}],'PixelRegion',{[1000 2000] [1000 2000]}));
    I = double(imread([SPN foldNames{i} tileNames{i}]));

    if ~exist('allI','var')
        allI = I * 0;
    end
    allI = allI + I/numTiles;
    %image(allI*(numTiles/i)), pause(.01)
    
end
save([SPN 'allI.mat'],'allI','-v7.3')
imwrite(uint8(allI),[SPN 'allI.tif'],'Compression','none')
        
    


               
            
            

