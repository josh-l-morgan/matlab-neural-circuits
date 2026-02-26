[TFN TPN] = GetMyFile('.mat')
nam = TFN(1:end-4);
load([TPN TFN])

secDir = [TPN nam];
if ~exist(secDir)
    mkdir(secDir);
end

for i = 1:size(sweeps,4)
    sprintf('Writing %d of %d.',i,size(sweeps,4))
    %image(uint8(sweeps(:,:,:,i)))
    imwrite(uint8(sweeps(:,:,:,i)),[secDir '\' num2str(i) '.tif']);
end