colormap(jet(64))
b = remez(10,[0 0.05 0.15 0.55 0.65 1],[0 0 1 1 0 0 ]);
[H,w] = freqz(b,1,128,'whole');
pause(1)
plot(2/pi-1,fftshift(abs(H)));
pause(1)


h = ftrans2(b);
freqz2(h)


bp1 = fir1(100,[.5 .80],'stop')
plot(bp1)

bp2 = fir2(48,[0 .1 .5 .9 1],[0 1 0 0 0]);
plot(bp2)
plot([0 .1 .5 .9 1],[0 1 0 0 0])

freak = [0 .1 .2 .3 .4 .5 .6 .7 .8 .8 1];
amp = [0 0 1 0 0 0 0 0 0 0 0];
myBP = fir2(50, freak, amp);
myBPr = ftrans2(myBP);
freqz2(myBPr)


%% 2D filter
fSize = 100;
kSize = 20;
for i = 1 :1
gKern = fspecial('gaussian',[kSize kSize],kSize/3);
gKern = gKern * 1/max(gKern(:));
gKern(fSize + kSize,fSize+kSize)=0;
image(fitH(gKern));


blah = rand(fSize,fSize);
blah = zeros(fSize,fSize);
blah(round(fSize/2),round(fSize/2)) = 1;
blah(1,1) = 1;
blah(fSize + kSize,fSize+kSize)=0;
fB = fft2(blah);
%image(fitH(fB))
fK = fft2(gKern);
%image(fitH(fK))

C = ifft2(fB.*fK);
C = ifft2(fft2(gKern).*fftn(blah))
C = real(C);
C = C(fix(kSize/2):fix(kSize/2)+fSize,fix(kSize/2):fix(kSize/2)+fSize);
image(fitH(C)),pause(.01)
[y x ] = find(C == max(C(:)))
end
