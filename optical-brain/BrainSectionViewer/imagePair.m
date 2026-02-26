function[Icol] = imagePair(p1, p2,ax);

m1 = mean(p1,3);
m2 = mean(p2,3);

[ys1 xs1] = size(m1);
[ys2 xs2] = size(m2);

ys = max(ys1,ys2);
xs = max(xs1,xs2);

Icol = zeros(ys,xs,3);

sy1 = fix((ys-ys1)/2);
sx1 = fix((xs-xs1)/2);
sy2 = fix((ys-ys2)/2);
sx2 = fix((xs-xs2)/2);

ms1 = m1-min(m1(:));
ms2 = m2-min(m2(:));

ms1 = ms1 * 50/median(m1(:));
ms2 = ms2 * 50/median(m2(:));

Icol(sy2+1:sy2+ys1,sx2+1:sx2+xs1,1) = ms2;
Icol(sy1+1:sy1+ys1,sx1+1:sx1+xs1,2) = ms1;
Icol(sy2+1:sy2+ys1,sx2+1:sx2+xs1,3) = ms2;

Icol = uint8(Icol);

image(ax,Icol)
xlim(ax,[1 size(Icol,2)])
ylim(ax,[1 size(Icol,1)])
drawnow



