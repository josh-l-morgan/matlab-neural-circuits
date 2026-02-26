


MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm_matOut\export_14+04+12_matOut\'
load([MPN 'dsObj.mat'])


fsize = [2000 2000 2000]
field = zeros(fsize,'uint16');
 
for i = 1:length(dsObj)
i

indI = sub2ind(fsize,dsObj(i).y,dsObj(i).x,dsObj(i).z);
field(indI) = i;

end

max1 = squeeze(max(field,[],1));



modMax = mod(max1,250);
modMax(max1>0) = modMax(max1>0)+2;
image(modField)

% 
% max2 = squeeze(max(field,[],2));
% max3 = squeeze(max(field,[],3));
% 



cmap = colormap(hsv(256));
cmap(1:2,:) = 0;
colormap(cmap)
image(mod(max1,200))
image(max2)
image(max3)


%%
I = zeros(fsize(1),fsize(2),'uint8');
for i = 1:length(dsObj)
i

indI = sub2ind([fsize(1) fsize(2)],dsObj(i).z,dsObj(i).x);
I(indI) = 1000;
image(I)
pause
I = I * 0;

end