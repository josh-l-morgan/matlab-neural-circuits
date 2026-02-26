colormap gray(256)
[TFN TPN] = GetMyFile;

I = imread([TPN TFN]);

image(I)
subplot(2,2,1)
image(I(900:1000,750:850))
subplot(2,2,1)
image(I)

If = medfilt2(I,[5 3]);
If = imfilter(If,fspecial('disk',21));

subplot(2,2,2)
image(If(900:1000,750:850))
subplot(2,2,2)
image(If)

subplot(2,2,3)
for i = 1:10: 256
    image((If>i)* 1000),pause% (.01);
   [Il n]= bwlabeln(If);
   n
end



