

[TFN TPN] = GetMyFile;  %Get file location
vidObj = mmreader([TPN TFN]);  %Get movie information as object
nFrames = vidObj.NumberOfFrames 
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;

nam = TFN(1:end-4);

sampFrames = 1;
subplot(1,1,1)

chunkSize = 1000;
readFrames = ((1:fix(nFrames/chunkSize))-1)*chunkSize;

%% Read avi

if exist([TPN nam '.tif'])
    delete([TPN nam '.tif'])
end

for i = readFrames
    tic
    sprintf('Reading frame %d of %d.', i,readFrames(end))
    frame = read(vidObj, [i+1 i+chunkSize]);  %read first few images
    for f = 1:size(frame,4)    
        %sprintf('Writing frame %d of %d.', f, size(frame,4))
        imwrite(frame(:,:,:,f),[TPN nam '.tif'],'writemode','append')
    end
    lastRead = i;
    chunckTime = toc;
    remainTime = chunkTime * sum(readFrames>i);
end

if nFrames>lastRead
    frame = read(vidObj, [lastRead+1 nFrames]);  %read first few images
    for f = 1:size(frame,4)    
        imwrite(frame(:,:,:,f),[TPN nam '.tif'],'writemode','append')
    end
end

