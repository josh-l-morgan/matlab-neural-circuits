clear all
[TFN DPN] = GetMyFile
TPN = [DPN TFN];

Iall = tiffread2(TPN);

zs = length(Iall);
xs = Iall(1).width;
ys = Iall(1).height;
planes = zs/2;
I = zeros(ys,xs,3,planes);

MedSize = 3;
for i = 1: zs

    I(:,:, fix((i-1)/planes)+1,mod(i-1,planes)+1) = medfilt2(Iall(i).data,[MedSize MedSize]);
end

%image(uint8(max(I,[],4)*2))

name = TFN(1:find(TFN =='.')-1)

for c = 1:2
I(:,:,c,:) = I(:,:,c,:) * 256/(max(max(max(I(:,:,c,:)))));
end

I = uint8(I);

Iwrite([DPN name],I)

'Done converting.'
