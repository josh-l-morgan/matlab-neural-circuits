%%Grab tiles again

%%
% SPN = 'G:\joshm\LGNs1\rawMontages\'
% TPN = 'E:\LGN_data\subVol2\grabedTiles\'
% rawDir = dir(SPN)

if ~exist(TPN,'dir'),mkdir(TPN),end

useWafer = {'w002','w003'}
%PixelRegion = {[10000 20000] [15000 25600]}

useRow = 2;
useCol = 2;


%%

pickSec = {}; c  = 0;
pickTile = {};
for w = 1:length(useWafer);
    currentWafer = useWafer{w};
    tilePre = sprintf('Tile_r%d-c%d_%s_',useRow,useCol,currentWafer);
    
    for i  = 1:length(rawDir)
        
        nam = rawDir(i).name;
        if rawDir(i).isdir
            if sum(regexp(nam,currentWafer))
                secDir = dir([SPN nam '\' tilePre '*.tif']);
                c = c+1;
                pickSec{c} = nam;
                pickTile{c} = secDir(1).name;
            end
        end
    end
end

%%

caughtError = 0;
for i = 1:length(pickSec)
    
    oldName = [SPN pickSec{i} '\' pickTile{i}];
    newName = [TPN pickTile{i}];
    
    if ~exist(newName,'file')
        disp(sprintf('Copying file %d of %d',i,length(pickSec)))
        
        
        if exist('PixelRegion')
            I = imread(oldName,'PixelRegion',PixelRegion);
            try
                imwrite(I,newName,'Compression','none')
            catch
                err
            end
        else
            try
                copyfile(oldName,newName)
            catch err
                caughtError = 1;
                err
            end
        end
    end
end

if caughtError
    disp('There was a copyfile error during transfer')
else
    disp('Copy completed without error.')
end



