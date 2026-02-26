clear all
colormap gray(256)

rad = 4;

[TFN TPN] = GetMyFile;

I = imread([TPN TFN]);
%I = 255 - I;
%I = I(1:1000,1:2000);
image(I),pause(.01)


for i = 1:1
   I = medfilt2(double(I),[3,3]); 
end
%%

[ys xs] = size(I);
vI = I;
cI = uint8(cat(3,I,(I)/3,I/3));

[y x] = find(I == min(I(:)),1);

lookShape = fspecial('disk',rad)>0;
[shifty shiftx] = find(lookShape);
shifty = shifty - rad - 1; shiftx = shiftx -rad - 1;
r = 0;
while 1
    r = r + 1;
    ny = shifty + y;
    nx = shiftx + x;
    ny(ny<1) = 1; ny(ny>ys) = ys;
    nx(nx<1) = 1; nx(nx>xs) = xs;
    
grab = sub2ind([ys xs], ny , nx);

gVals = vI(grab);
go2 = grab(find(gVals == min(gVals)));
goTo = go2(fix(rand * length(go2))+1);
vI(goTo) = vI(goTo) + 1;
vI(grab) = vI(grab) + 1;

[y x] = ind2sub([ys xs], goTo);
if ~(mod(r,300))
    tI = (vI-I) ;
    tI = tI * 50/median(tI(tI(:)>0));
    tI(tI>0) = tI(tI>0) + 50;
    tI(grab) = 1000;
    cI(:,:,1) = tI;
    
    image(cI),pause(.01)
end
end

cI2 = cI;
cI2(:,:,1) = 255 - cI2(:,:,3);
cI2(:,:,2) = 255 - cI2(:,:,2);
image(cI2)

save([TPN 'path.mat'],'cI')

image(((vI-I)>70)*1000)

image(255-(vI-I) * 2)
% lI = bwlabeln(