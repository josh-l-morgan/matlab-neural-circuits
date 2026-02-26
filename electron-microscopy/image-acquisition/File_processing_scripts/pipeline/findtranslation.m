function [transx,transy,transcorr]=findtranslation(image1,image2,filtmin,filtmax,baseimg)
%A function to compute the best-fit translation between images 1 and 2
%filtmin and filtmax give the minimal and maximal periods of spatial frequencies
%which should remain in the image (in pixels).
%baseimg gives the image number from which processing images should be rendered;
%give 0 here to suppress image output.
%By Daniel Berger, March 18, 2009 (for MIT-BCS Seung)

fim1=bandpass2(image1,filtmin,filtmax);
if baseimg>0
  figure(baseimg);
  vim=(fim1-min(min(fim1)))/(max(max(fim1))-min(min(fim1)));
  imshow(vim);
  title('Filtered cut image 1');
end;

fim2=bandpass2(image2,filtmin,filtmax);
if baseimg>0
  figure(baseimg+1);
  vim=(fim2-min(min(fim2)))/(max(max(fim2))-min(min(fim2)));
  imshow(vim);
  title('Filtered cut image 2');
end;

%transmap=xcorr2(fim1,fim2);
transmap=convn_fast(fim1,flipdims(fim2),'full');
if baseimg>0
  figure(baseimg+2);
  imagesc(transmap);
  title('convn_fast translation map');
end;
            
%compute displacement vector
middx=floor(size(transmap,2)/2+1);
middy=floor(size(transmap,1)/2+1);
[A,B]=max(transmap);
[C,D]=max(A);
transx=D-middx;  %X-translation of second image leading to maximal correlation
transy=B(D)-middy; %y-translation of second image leading to maximal correlation
transcorr=C;  %value of the maximal correlation

if baseimg>0
  figure(baseimg+2);
  hold on;
  plot(D,B(D),'wo');
  hold off;
end;
