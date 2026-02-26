function[kern] = addKern(kern1,kern2);


[ys1 xs1] = size(kern1);
[ys2 xs2] = size(kern2);

ys3 = max(ys1,ys2);
xs3 = max(xs1,xs2);

kern = zeros(ys3,xs3);

kern((ys3-ys1)/2+1:ys3-(ys3-ys1)/2,(xs3-xs1)/2+1:xs3-(xs3-xs1)/2) = ...
    kern((ys3-ys1)/2+1:ys3-(ys3-ys1)/2,(xs3-xs1)/2+1:xs3-(xs3-xs1)/2) + kern1;

kern((ys3-ys2)/2+1:ys3-(ys3-ys2)/2,(xs3-xs2)/2+1:xs3-(xs3-xs2)/2) = ...
    kern((ys3-ys2)/2+1:ys3-(ys3-ys2)/2,(xs3-xs2)/2+1:xs3-(xs3-xs2)/2) + kern2;

%image(kern*100)