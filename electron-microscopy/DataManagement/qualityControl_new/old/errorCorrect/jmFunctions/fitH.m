function[Im] = fitH(I,p)

iMin = double(min(I(:)));
iMax = double(max(I(:)));


I = (I - iMin)* 255/(iMax-iMin);

if nargin==2
    Im = I(:,:,p);
else
    Im = I;
end