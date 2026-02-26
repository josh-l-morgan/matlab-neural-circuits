

SPN = 'G:\grant_images\hxQ\60x_target\MatProc\2\'
TPN = 'G:\grant_images\hxQ\60x_target\MatProc\2_gaus\'
colormap gray(256)
if ~exist(TPN,'dir'), mkdir(TPN),end
Inam = dir([SPN '*.tif']);

It = imread([SPN Inam(1).name]);

[ys xs cs] = size(It);

look = 2;

Iv = zeros(ys,xs,length(Inam),'double');
for i = 1:length(Inam)
    Iv(:,:,i) = imread([SPN Inam(i).name]);
end


if 0 % med 3d
    
    Im = medfilt3(Iv,[5 5 9]);
    
else
    
    Im = imgaussfilt3(Iv*20,6);
end

%Im = Im * 255/max(Im(:));
for i = 1:size(Im,3)
    image(Im(:,:,i)),pause(.1)
    imwrite(uint8(Im(:,:,i)),sprintf('%smed_%04.0f.tif',TPN,i))
end
%
% for i = 1:length(Inam)
%    range = [max(1,i-look): min(length(Inam),i+look)];
%
%    Ibuf = zeros(ys, xs, length(range),'double');
%    
%    for r = 1:length(range);
%        Ibuf(:,:,r) = imread([SPN Inam(range(r)).name]);
%    end
%    
%    Im = median(Ibuf,3);
%     
%    subplot(1,2,1)
%     image(Im), 
%     If = medfilt
%     subplot(1,2,2)
%     image(Im), 
% 
%     pause(.1)
%     %image(Ibuf(:,:,1))
% end


































