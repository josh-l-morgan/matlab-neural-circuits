function[kern] = cone(kSize);




kern = ones(kSize);
[y x] = find(kern);
dists = sqrt((y-kSize/2).^2 + (x-kSize/2).^2);
kern(:) = dists;
cutVal = kern(round(kSize/2),1);
kern = cutVal - kern;
kern(kern<0) = 0;

image(kern*256/max(kern(:)));
