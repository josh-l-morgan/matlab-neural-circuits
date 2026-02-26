function [c,i]=getrotation(image, startang, step, endang, hpfcutoff, basefig)
%getrotation.m

%close all;
%clear all;
%image=imread('slice1024_01.png');
%image=imread('hip29nm_W4S11R1C1_10x.png');

%image=image(:,:,1);

%Show image
if basefig>0
  figure(basefig);
  imshow(image);
  title('Input image');
end;

%figure(2);
%im2=fft2(image);
%surf(squeeze(log(abs(im2(1:10:size(im2,1),1:10:size(im2,2),3)))));

%Compute and show Radon transform
theta = startang:step:endang;
[R,xp] = radon(image,theta);

if basefig>0
  figure(1+basefig);
  imagesc(theta,xp,R);
  title('Radon transform of input image');
end;

%Compute and show Fourier transform of Radon transform
im=fft(R);

if basefig>0
  figure(2+basefig);
  imagesc(log(abs(im)));
  title('Fourier transform of Radon transform');
end;

%Compute and show subregion R2 of radon transform
%b=floor(size(R,1)/2-min(size(image))/2);
%e=floor(size(R,1)/2+min(size(image))/2);
b=floor(size(R,1)/2-min(size(image))/4);
e=floor(size(R,1)/2+min(size(image))/4);
R2=R(b:e,:);

if basefig>0
  figure(3+basefig);
  imagesc(R2);
  title('Subimage of Radon transform Image, without edges');
end;

im2=fft(R2);

if basefig>0
  figure(4+basefig);
  %imagesc(log(abs(im2)));
  imagesc(log(abs(im2(hpfcutoff:1:size(im2,1)-hpfcutoff,:))));
  title('Fourier transform of subimage of Radon transform');
end;

%mim2=max(abs(im2(10:1:size(im2,1)-10,:)));
%[c,i]=max(mim2)

mim3=sum(abs(im2(hpfcutoff:1:size(im2,1)-hpfcutoff,:)));

if basefig>0
  figure(5+basefig);
  plot(theta,mim3);
  grid on;
  title('Sum of higher frequencies of Fourier spectrum (estimated slice orientation)');
end;

[c,i]=max(mim3);
i=i*step+startang;