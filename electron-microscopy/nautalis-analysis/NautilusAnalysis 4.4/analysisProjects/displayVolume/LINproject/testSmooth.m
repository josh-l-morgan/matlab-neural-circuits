

%%
colormap gray(256)
f = zeros(101,101,101);
f(51,51,51) = 100;

k = 51;
s = smooth3(f,'gaussian',k,k/5);


subplot(2,1,1)
sumS = squeeze(sum(s,3));
image(sumS * 256/max(sumS(:)));
subplot(2,1,2)
plot(s(51,:,51))