

fileName = 'F:\cuttingVideos\20151119_110913_testPSB_kn45_11-18_startCortex_40nm_IonizerBigClose.avi'
xyloObj = VideoReader(fileName);
info = get(xyloObj)



nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

       
 mov1 = zeros(vidHeight,vidWidth,nFrames,3);      
% Read one frame at a time.
for k = 1 : nFrames
    mov(k).cdata = read(xyloObj, k);
    mov1(:,:,i,:) = mov(i).cdata;
end

% Play back the movie once at the video's frame rate.
movie(mov, 1, xyloObj.FrameRate);

for i = 1:length(mov)
    mov1(:,:,i,:) = mov(i).cdata;
end

for i = 1 : size(mov1,3)
   image(mov1(:,:,i,:)); 
end