function[kern] = tent(kSize,fallOff);


%%
if ~exist('fallOff')
    fallOff = 2;
end

if ~mod(kSize,2)
    kSize = kSize+1;
    'making kernal odd'
end


kern = ones(kSize);
[y x] = find(kern);
dists = sqrt((y-(kSize+1)/2).^2 + (x-(kSize+1)/2).^2);
kern(:) = dists.^(1/fallOff);
cutVal = kern(round(kSize/2),1);
kern = cutVal - kern;
kern(kern<0) = 0;


%{
    [Gx,Gy] = imgradientxy(kern,'sobel');


subplot(2,1,1)
bar(kern(round(kSize/2),:))
subplot(2,1,2)
bar(Gx(round(kSize/2),:))

d = diag(Gx)
bar(d)
%}
