

colormap gray(256)
TPN = GetMyDir;

nams = getPics(TPN);
info = imfinfo([TPN nams{1}]);

for i = 1: length(nams)
    I = double(imread([TPN nams{i}]));
    gKern2 = gaus3d([200 200 1],1); 
    gKern1 = gaus3d([200 200 1],10); %df = ([ 20 20 1], 3)
    gKern3 = gKern2-gKern1;
    image(gKern3)
    c1 = fastCon2d(I,gKern3);
    image(c1)
    %fit
    c1(c1>mean(c1(:))) = mean(c1(:));
    c1 = c1-min(c1(:));
    c1 = c1/(max(c1(:)))*255;
    image(c1)
    
%     
% h = fspecial('laplacian', 20)
% J = real(ifft2(fft2(I).*fft2(h, size(I, 1), size(I, 2))));
% image(J)
%     
        
    
    
    
end



