


kern1 = ones(3, 3, 3);
rad = 30;
testVol = zeros(rad*2,rad*6,rad*2);
testVol(rad,rad,rad) = 1;
If = convn(testVol,kern1);



kern1 = ones(7, 7, 7);
rad = 30;
testVol = zeros(rad*2,rad*6,rad*2);
testVol(1:rad*2,rad,rad) = 1;
testVol(rad,rad,1:rad*2) = 1;

If = convn(testVol,kern1);

colormap gray(256)
sumI = sum(double(If),3);
image(sumI * 256/max(sumI(:)))

indI = find(If);
[y x z] = ind2sub(size(If),indI);
subs = [y x z];
seed = subs(1,:);

subs2arbor(subs,seed)