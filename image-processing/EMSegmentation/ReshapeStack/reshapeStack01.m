


SPN = 'Z:\Active\morganLab\DATA\LGN\LGNs1\tweakedImageVolume\LGNs1_P32_smallHighres\'
TPN = 'Z:\Active\morganLab\DATA\LGN\LGNs1\tweakedImageVolume\TestDoubling\'

SPNd = dir([SPN 'Segmentation1*.png']);
sNams = {SPNd.name};

SPNd = dir([SPN 'tweakedImageVolume2*.png']);
iNams = {SPNd.name};


bDepth = 256;
hsvMap = hsv(1000);
cmap = hsvMap(randperm(1000),:);




for s = 1:length(sNams)

    snam = sNams{s};
    Is = double(imread([SPN snam]));
    Isi = Is(:,:,3) + Is(:,:,2) * bDepth + Is(:,:,1) * bDepth ^2;
    IsiD = imresize(Isi,.5,'nearest');
    IsD = imresize3(Is,[size(Is,1)/2 size(Is,2) size(Is,3)],'nearest');

    inam = iNams{s};
    I = imread([SPN inam]);

    z = (s-1) * 2 + 1;
    z2 = (s-1) * 2 + 2;

    imwrite(uint8(IsD),[TPN sprintf('seg%04.0f.png',z)]);
    imwrite(I,[TPN sprintf('sec%04.0f.png',z)]);
    imwrite(uint8(IsD),[TPN sprintf('seg%04.0f.png',z2)]);
    imwrite(I,[TPN sprintf('sec%04.0f.png',z2)]);


end







