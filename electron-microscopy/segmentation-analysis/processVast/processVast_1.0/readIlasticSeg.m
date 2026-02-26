
SPN = 'D:\LGNs1\segmentation\Odyssey2_trial\32nm_40-41\Ilastik\image2\segmentation.png'


I = double(imread(SPN));
cyto = I(:,:,1)./I(:,:,2)>1;

labI = bwlabel(cyto);

image(mod(labI,100))