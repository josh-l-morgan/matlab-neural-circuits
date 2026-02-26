function[I2] = outline(I,r);

if ~exist('r','var')
    r = 5;
end

bg = 0;
backG = I==bg;
mask = sum(backG,3) < 3;

se = strel('disk',r);

di = imdilate(mask,se);

I2 = I;

I2(repmat((di>0) & ~mask,[1 1 3])) = .001;
% 
% image(If)
% image(I)
