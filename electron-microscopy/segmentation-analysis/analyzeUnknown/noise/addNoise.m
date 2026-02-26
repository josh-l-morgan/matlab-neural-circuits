
colormap gray(256)
[IFN IPN] = GetMyFile
[NFN NPN] = GetMyFile
[FFN FPN] = GetMyFile

I = imread([IPN IFN]);
N = imread([NPN NFN]);
F = imread([FPN FFN]);


I = double(I);
N = double(N);
P = I + N-mean(N(:));

subplot(2,2,1)
image(I)
subplot(2,2,2)
image(N)
subplot(2,2,3)
image(P)
subplot(2,2,4)
image(F)