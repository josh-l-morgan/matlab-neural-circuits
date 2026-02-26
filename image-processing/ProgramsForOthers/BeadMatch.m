[TFN TPN] = GetMyFile;
%tiffread2([TPN TFN]);

fname = [TPN TFN];
info = imfinfo(fname);
num_images = numel(info);
for k = 1:num_images
    A = imread(fname, k);
    I(:,:,k) = double(A(:,:,1));
end

Icon = I;
%Icon = I(450:500,150:250,:);
gkern = gaus3d([10,10,3],3);
Icon = fastCon(Icon,gkern);
subplot(1,2,1)
image(max(I,[],3));
subplot(1,2,2)
image(max(Icon,[],3)*10);
pause(.01)

%%
wI = max(Icon(:))-Icon;
wI = imhmin(wI,5);
wI(wI>=mean(wI(:))) = 0;
wI = watershed(wI,26);

colormap colorcube
for i = 1:size(I,3)
image(wI(:,:,i)), pause
end

%%
mapI = wI * 0;
for i = 2:max(wI(:))
   inds = find(wI==i);
   [y x z] = ind2sub(size(wI),inds); 
   vals = Icon(inds);
   b(i).meanBright = mean(vals);
   b(i).size = length(vals);
   vals = vals-mean(vals);
   vals(vals<0) = 0;
   vals = vals.^3;
   tot = sum(vals);
   my = sum((y .* vals))/tot;
   mx = sum((x .* vals))/tot;
   mz = sum((z .* vals))/tot;
   b(i).Y = my;
   b(i).X = mx;
   b(i).Z = mz;
   mapI(inds) = 100;
   mapI(round(my),round(mx),round(mz)) = 1000;
   
end


