colormap gray(256)
[TFN TPN] = GetMyFile;

I = imread([TPN TFN]);
I = I * 256/max(I(:));
Ir = double(I(:,:,1));
Ig = double(I(:,:,2));
subplot(2,2,1)
image(255 - Ir)
subplot(2,2,2)
image(255 - Ig)


[ys, xs] = size(Ir);

%% auto Correlation
sRange = 10;
clear cc1 cc2
for y = -sRange :1: sRange
    for x = -sRange :1 : sRange
        grabIg = Ig( sRange + 1 + y:ys - sRange + y, sRange + 1 + x:xs - sRange + x) ;
        
        cc = corrcoef(grabIg,Ig(sRange + 1:ys - sRange,sRange + 1:xs - sRange));
        cc2(sRange + 1 + y,sRange + 1 +x) = cc(2,1);
        
        grabIr = Ir( sRange + 1 + y:ys - sRange + y, sRange + 1 + x:xs - sRange + x) ;
        
        cc = corrcoef(grabIr,Ir(sRange + 1:ys - sRange,sRange + 1:xs - sRange));
        cc1(sRange + 1 + y,sRange + 1 +x) = cc(2,1);
    end
end
subplot(2,2,3)
image(cc1 * 255/max(cc1(:)))
subplot(2,2,3)
image(cc2 * 255/max(cc2(:)))

corDif = cc2./cc1;
%image(corDif * 255/max(corDif(:))),pause(.01)

xfilt = cc2(sRange+ 1,:);
%%
% 
for l = 2:7
 [J , PSF] = deconvblind(Ig,ones(1,l));
% J = deconvreg(Ig, xfilt,.6);
% J = J - min(J(:));
% J = J * 255/max(J(:));
subplot(2,2,4)
image(255 - J),pause
% 
% subplot(2,2,3)
% image(255 - J) 
end
%%


fminsearch

%imwrite(uint8(255-If),[TPN 'Filt' TFN],'tif','Compression','none')