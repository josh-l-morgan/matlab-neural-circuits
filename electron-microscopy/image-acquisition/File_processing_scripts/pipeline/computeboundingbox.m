function bounds=computeboundingbox(imgx1,imgy1,imgx2,imgy2, tmatrix)
%A function to compute the bounding box of an image with given 
%(imgx1,imgy1)-(imgx2,imgy2) position and dimensions after it is transformed 
%by the given 3x3 tmatrix.
%Returns four values in the following order: [minx, miny, maxx, maxy]
%By Daniel Berger for MIT-BCS Seung, June 1st 2009

xy1=tmatrix*[imgx1 imgy1 1]';
xy2=tmatrix*[imgx2 imgy1 1]';
xy3=tmatrix*[imgx2 imgy2 1]';
xy4=tmatrix*[imgx1 imgy2 1]';

bounds=zeros(1,4);
bounds(1)=min([xy1(1) xy2(1) xy3(1) xy4(1)]);
bounds(2)=min([xy1(2) xy2(2) xy3(2) xy4(2)]);
bounds(3)=max([xy1(1) xy2(1) xy3(1) xy4(1)]);
bounds(4)=max([xy1(2) xy2(2) xy3(2) xy4(2)]);