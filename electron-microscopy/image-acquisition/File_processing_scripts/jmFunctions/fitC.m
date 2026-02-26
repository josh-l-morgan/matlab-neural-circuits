function[Im] = fitC(rawI,p)

Im = rawI;
for i = 1:size(rawI,3)
    iMin = double(min(I(:)));
    iMax = double(max(I(:)));


    I = (I - iMin)* 255/(iMax-iMin);


    Im(:,:,i) = I;
end
colormap gray(255)