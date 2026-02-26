function bounds=computeboundingbox_recentered(imgx1,imgy1,imgx2,imgy2, tmatrix)
%A function to compute the bounding box of an image with given 
%(imgx1,imgy1)-(imgx2,imgy2) position and dimensions after it is transformed 
%by the given 3x3 tmatrix.
%This version uses re-centering placing the affine transform at the center
%of the image.
%Returns four values in the following order: [minx, miny, maxx, maxy]
%By Daniel Berger for MIT-BCS Seung, June 1st 2009

midx=(imgx2+imgx1)/2;
midy=(imgy2+imgy1)/2;

xy1=tmatrix*[imgx1-midx imgy1-midy 1]';
xy2=tmatrix*[imgx2-midx imgy1-midy 1]';
xy3=tmatrix*[imgx2-midx imgy2-midy 1]';
xy4=tmatrix*[imgx1-midx imgy2-midy 1]';

bounds=zeros(1,4);
bounds(1)=min([xy1(1) xy2(1) xy3(1) xy4(1)])+midx;
bounds(2)=min([xy1(2) xy2(2) xy3(2) xy4(2)])+midy;
bounds(3)=max([xy1(1) xy2(1) xy3(1) xy4(1)])+midx;
bounds(4)=max([xy1(2) xy2(2) xy3(2) xy4(2)])+midy;