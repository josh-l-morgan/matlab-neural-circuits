%%subtract template
colormap gray(256)
OPN = GetMyDir;
nams = getPics(OPN);

TPN = GetMyDir;
tnams = getPics(TPN);

I = imread([TPN tnams{1}]);
Isum = double(medfilt2(I,[3,3]));
for i = 2:length(tnams)
    I = imread([TPN tnams{i}]);
    I = medfilt2(I,[3,3]);
    Isum = Isum + double(I);
    i
end
Imean = Isum/length(tnams);
image(I)


if ~exist([OPN(1:end-1) 'fix\'])
    mkdir([OPN(1:end-1) 'fix\'])
end

for i = 1:length(nams)
    I = imread([OPN nams{i}]);
    Inew = I + uint8(fix);
    subplot(1,2,1)
    image((I))
    subplot(1,2,2)
    image((Inew))
  %imwrite(Inew,[TPN(1:end-1) 'fix\' num2str(i) '.tif'],'Compression','none')
  imwrite(Inew,[OPN(1:end-1) 'fix\' char(nams{i})],'Compression','none')
   i
end







