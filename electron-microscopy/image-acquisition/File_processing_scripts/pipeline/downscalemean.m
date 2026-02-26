function rimg=downscalemean(img, blocksize)
%takes a gray image img and downsamples it to rimg by using the min in each
%block. Image size has to be a divisible by blocksize.

rimg=zeros(size(img,1)/blocksize,size(img,2)/blocksize);

ty=1;
for y=1:blocksize:size(img,1)
  tx=1;
  for x=1:blocksize:size(img,2)
    rimg(ty,tx)=mean(mean(img(y:y+blocksize-1,x:x+blocksize-1)));
    tx=tx+1;
  end;
  ty=ty+1;
end;