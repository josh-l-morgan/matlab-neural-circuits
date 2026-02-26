
clear all
fileName = 'C:\Users\joshm\Documents\myWork\myPublications\LGNs1\figures\movies\rotation\combineSeedAssociation\cropI1small.avi'
xyloObj = VideoReader(fileName);
info = get(xyloObj)



nFrames = xyloObj.NumberOfFrames;
vidHeight = xyloObj.Height;
vidWidth = xyloObj.Width;

% Preallocate movie structure.
mov(1:nFrames) = ...
    struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
           'colormap', []);

       
 mov1 = zeros(vidHeight,vidWidth,3,nFrames,'uint8');      
% Read one frame at a time.
parfor k = 1 : nFrames
    mov(k).cdata = read(xyloObj, k);
end

% Play back the movie once at the video's frame rate.
%movie(mov, 1, xyloObj.FrameRate);

SE = strel('disk',1);
clear mov1
parfor i = 1:length(mov)
    If = double(mov(i).cdata);
     gamma = .6;
    If = If.^gamma * 256/256^gamma;
    %If = imdilate(If,SE);
    If = imresize(If,.75);
    mov1(:,:,:,i) = uint8(If);
end

% for i = 1 : size(mov1,4)
%    image(squeeze(mov1(:,:,:,i))); 
%    pause(.01)
% end

dotAv = regexp(fileName,'.avi');
compName =[fileName(1:dotAv-1) 'Compress8.avi'];

myVideo = VideoWriter(compName);
myVideo.FrameRate = 30;  % Default 30
myVideo.Quality = 25;    % Default 75
open(myVideo);
writeVideo(myVideo, uint8(mov1));
close(myVideo);



