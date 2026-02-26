clear all
colormap gray(256)
[TFN TPN] = GetMyFile;



I = double(imread([TPN TFN]));
I = I * 256/max(I(:));
image(I)

[J , PSF] = deconvblind(I,ones(1,7));

image(J)

%%  Align


%imwrite(uint8(255-If),[TPN 'Filt' TFN],'tif','Compression','none')