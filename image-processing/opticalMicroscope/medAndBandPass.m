

SPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack_tweak6\'
TPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack_tweak6Gau\'
colormap gray(256)

useDif = 0;


if ~exist(TPN,'dir'), mkdir(TPN),end
Inam = dir([SPN '*.tif']);

It = imread([SPN Inam(1).name]);

It = imresize(It,1/4);
[ys xs cs] = size(It);

look = 2;

zs = length(Inam);
Ir = zeros(ys,xs,length(Inam),'double');
Ig = Ir;
Ib = Ir;
for i = 1:length(Inam)
    I = double(imread([SPN Inam(i).name]));
    I = imresize(I,1/4);
    Ir(:,:,i) = I(:,:,1);
    Ig(:,:,i) = I(:,:,2);
    Ib(:,:,i) = I(:,:,3);
end

if useDif
    Irf = imgaussfilt3(Ir,1);
    Irf2 = imgaussfilt3(Ir,16);
    Ird = Irf2-Irf;
else
    Ird = imgaussfilt3(Ir,1);
end


for i = 1:zs
    image(Ird(:,:,i)*1255/max(Ird(:)))
    pause(.3)
end

if useDif
    
    Igf = imgaussfilt3(Ig,1);
    Igf2 = imgaussfilt3(Ig,16);
    Igd = Igf2-Ig;
else
    Igd = imgaussfilt3(Ig,1);
end
for i = 1:zs
    image(Igd(:,:,i)*255/max(Igd(:)))
    pause(.01)
end

if useDif
    
    Ibf = imgaussfilt3(Ib,1);
    Ibf2 = imgaussfilt3(Ib,16);
    Ibd = Ibf2-Ibf;
else
    Igd = imgaussfilt3(Ig,1);
end
for i = 1:zs
    image(Ibd(:,:,i)*255/max(Ibd(:)))
    pause(.01)
end




%Im = Im * 255/max(Im(:));
for i = 1:zs
    Ic = cat(3,Ird(:,:,i),Igd(:,:,i),Ibd(:,:,i));
    Ic = Ic * 2^8/2^16;
    imwrite(uint8(Ic),sprintf('%sbp_%04.0f.tif',TPN,i))
end
































