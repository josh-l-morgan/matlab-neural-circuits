%%Grab tiles again

%%
SPN = 'E:\LGN_data\subVol2\grabedTiles\'
TPN = 'E:\LGN_data\subVol2\grabedTilesSub\'

rawDir = dir(SPN)

if ~exist(TPN,'dir'),mkdir(TPN),end

PixelRegion = {[10080 19328] [15168 24672]}

%%
c  = 0;
pickTile = {};

for i  = 1:length(rawDir)
    pickTile{i} =rawDir(i).name;
end

%%

caughtError = 0;
for i = 1:length(pickTile)
    
    oldName = [SPN  pickTile{i}];
    newName = [TPN pickTile{i}];
    
    if ~exist(newName,'file')
        disp(sprintf('Copying file %d of %d',i,length(pickTile)))
        
        
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



