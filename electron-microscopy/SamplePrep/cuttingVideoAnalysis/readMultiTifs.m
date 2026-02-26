

[TFN TPN] = GetMyFile


Info = imfinfo([TPN TFN]);
fSize = length(Info);

for i = 1:fSize
   lscan = imread([TPN TFN],'Index',i,'PixelRegion',{[2 200],[4 300]})
end
