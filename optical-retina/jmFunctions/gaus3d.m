function[fKern] = gaus3d(kSize,kSig,pixSize)

%%Make a 3d gaussian kernal with the dimensions kSize
%%the sigma of kSig 
%%the size of the pixels can be specified with pixSize, otherwise asumes 
%%1 x 1 x 1

if nargin <3
    pixSize = [1 1 1];
end

kRad = (kSize + 1)/2;
fKern = zeros(kSize);

[y x z] = ind2sub(kSize,find(fKern==0));
dists = sqrt(((y-kRad(1))*pixSize(1)).^2 ...
    + ((x - kRad(2))*pixSize(2)).^2 + ((z-kRad(3))*pixSize(3)).^2);

fKern(:) = 1 * exp(-.5 * (dists/kSig).^2);
fKern = fKern/sum(fKern(:));
