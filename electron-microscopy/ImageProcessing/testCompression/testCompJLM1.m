
colormap gray(256)
fileName = 'C:\Users\joshm\Documents\myWork\myProjects\Alignment\HPrealign\Compression\Tile_r2-c1_w041_sec016-1_NiceImage.tif';


Iraw = imread(fileName);



useWin = [3401 3600; 5401 5600];
useWin = [3001 3600; 5001 5600];

I = Iraw(useWin(1,1):useWin(1,2),useWin(2,1):useWin(2,2));
image(I)

%%
medKern = [1 1 1; 1 1 1; 1 1 1];
[yshift xshift] = find(medKern);
yshift = yshift-2;
xshift = xshift-2;

vertKern = 1:6;
twoKern = [1 2 4 5];

If1 = I;
If2 = I;
If3 = I;
for y = 2:size(I,1)-1
    sprintf('running line %d of %d',y,size(I,1)-1)
    parfor x = 2:size(I,2)-1
        vind = sub2ind(size(I),y+yshift,x+xshift);
        val = double(I(vind));
        If1(y,x) = median(val);
        If2(y,x) = mean(val);
        val2 = double(I(vind(twoKern)));
        If3(y,x) = median(val2);
        
    end
end


%%

subplot(2,2,1) 
image(I)
title('raw')

subplot(2,2,2) 
image(If1)
title('3 x 3 median ')
subplot(2,2,3)
image(If2)
title('3 x 3 mean ')

subplot(2,2,4)
image(If3)
title('3 x 2 median ')


%% 
saveFold = 'C:\Users\joshm\Documents\myWork\myProjects\Alignment\HPrealign\Compression\compTest\';
%mkdir(saveFold);
saveQual = 30;
imwrite(uint8(I),[saveFold 'I.jpg'],'Quality',saveQual)
imwrite(uint8(If1),[saveFold 'If1.jpg'],'Quality',saveQual)
imwrite(uint8(If2),[saveFold 'If2.jpg'],'Quality',saveQual)
imwrite(uint8(If3),[saveFold 'If3.jpg'],'Quality',saveQual)

Ij = imread([saveFold,'I.jpg']);
If1j = imread([saveFold,'If1.jpg']);
If2j = imread([saveFold,'If2.jpg']);
If3j = imread([saveFold,'If3.jpg']);




subplot(2,2,1) 
image(Ij)
title('raw')

subplot(2,2,2) 
image(If1j)
title('3 x 3 median ')
subplot(2,2,3)
image(If2j)
title('3 x 3 mean ')

subplot(2,2,4)
image(If3j)
title('3 x 2 median ')



%% Fast FFt

kern  = [0 1 0; 1 2 1; 0 1 0];
kern = kern/sum(kern(:));

Ic1 = fastCon(I,kern);




