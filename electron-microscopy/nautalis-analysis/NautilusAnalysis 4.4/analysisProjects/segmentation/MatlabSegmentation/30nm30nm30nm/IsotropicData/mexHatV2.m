function[Imf] = mexHatV2(I,s1,s2);

%% process I
kSize = [max(s1,s2)*3 max(s1,s2)*3];
kRad = (kSize + 1)/2;
kern = zeros(kSize);

[y x z] = ind2sub(kSize,find(kern==0));
dists = sqrt(((y-kRad(1))).^2 + ((x - kRad(2))).^2);

cKern = 1 * exp(-.5 * (dists/s1).^2);
cKern = cKern/sum(cKern(:));
sKern = 1 * exp(-.5 * (dists/s2).^2);
sKern = sKern/sum(sKern(:));

if s2
kern(:) = cKern - sKern;
else
   kern(:) = cKern ; 
end

% subplot(2,1,1)
% plot(kern(round(kRad(1)),:))

%% Convolve

Itemp = fastCon(I,kern);

pixClip = s2;
Imf  = Itemp * 0;
Imf(pixClip +1:end-pixClip,pixClip +1:end-pixClip)=Itemp(pixClip +1:end-pixClip,pixClip +1:end-pixClip);
Imf(Imf<0)=0;
