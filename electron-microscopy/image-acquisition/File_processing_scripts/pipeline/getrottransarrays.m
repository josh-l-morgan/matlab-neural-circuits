function [transx,transy,corrval]=getrottransarrays(img1,img2,angles,method,baseimage)
%If baseimage is 0, no images will be shown.
%Otherwise, images starting at baseimage will be shown.

%angstep=(endangle-startangle)/nrofangles*pi/180;
nrofangles=max(size(angles));

%rots=[1:1:nrofangles]*(endangle-startangle)/nrofangles;
corrval=zeros(1,nrofangles);
transx=zeros(1,nrofangles);
transy=zeros(1,nrofangles);

s=floor(size(img2,1)*sqrt(2)/2-2);
for rot=1:1:nrofangles
  %rot
  rimg2=rotateimage(img2,angles(rot)*pi/180,s,s,method);
  rimg2=rimg2-mean(mean(rimg2)); %re-normalize
  
  transmap=convn_fast(img1,flipdims(rimg2),'full');
  %transmap=xcorr2(img1,rimg2);
  
  if baseimage>0
    figure(baseimage);
    imagesc(img1);
    colormap(gray);
    axis square;
    title('Image 1');
    figure(baseimage+1);
    imagesc(rimg2);
    colormap(gray);
    title('Rotated image 2');
    axis square;
    figure(baseimage+2);
    imagesc(transmap);
    title('Cross-correlation for different translations');
  end;

  middx=floor(size(transmap,2)/2+1);
  middy=floor(size(transmap,1)/2+1);
  [A,B]=max(transmap);
  [C,D]=max(A);
  transx(rot)=D-middx;  %X-translation of second image leading to maximal correlation
  transy(rot)=B(D)-middy; %y-translation of second image leading to maximal correlation
  corrval(rot)=C;
  %pause(0.1);
  %input('Press return...');
end;