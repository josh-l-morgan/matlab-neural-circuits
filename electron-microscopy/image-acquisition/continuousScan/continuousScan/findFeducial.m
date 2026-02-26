% 
% colormap gray(256)
% [IFN IPN] = GetMyFile;
% [FFN FPN] = GetMyFile;
% I = imread([IPN IFN]);
% F = imread([FPN FFN]);
% subplot(2,1,1)
% image(I)
% subplot(2,1,2)
% image(F)

%function[] = findFeducial(I,F)

%% 
sI = single(imresize(I,.1));
subplot(2,1,1)
image(fitH(sI))
sF = single(imresize(F,.1));
subplot(2,1,2)
image(sF)

%% filter sI
dSize = 3;
aKern = fspecial('disk',dSize);
aKern(dSize+1,dSize+1) = 0;
aKern = aKern/sum(aKern(:));
aKern(dSize+1,dSize+1) = -1;
sum(aKern(:))

dI = fastCon(sI,aKern);
subplot(2,1,1)
image(dI+100)

%% filter sF

dF = fastCon(sF,aKern);
subplot(2,1,2)
image(abs(dF))

%%
kF = dF(dSize:end-dSize,dSize:end-dSize);
fc = fastCon(abs(dI),abs(kF));
image(fitH(fc))

%% Get max
fI = fc== max(fc(:));
fI = fastCon(fI,fspecial('disk',10));
image( fitH(fI) )



