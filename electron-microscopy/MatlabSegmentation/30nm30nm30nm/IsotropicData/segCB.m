

colormap(gray(256))
SPN = GetMyDir;
TPN = [SPN(1:end-1) '_colorCBV' '\'];
mkdir(TPN)
%% read names
dSPN = dir(SPN);
iNam = {};
for i  =   1:length(dSPN)
    nam = dSPN(i).name;
    if sum(regexp(nam,'.tif'))
        iNam{length(iNam)+1} = nam;
    end
end

%%


colI = zeros(800,800,length(iNam),3,'uint8');

for i = 1:length(iNam)
   colI(:,:,i,:) = imread([SPN iNam{i}]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
%%Run for cell bodies

I = colI(:,:,:,2)>0;
G = colI(:,:,:,3);

%% 
[labI lnum] = bwlabeln(I,6);
iProps = regionprops(labI,'PixelIdxList');




%% 
comI3 = labI*0;
comI2 = comI3;
for o = 1:lnum
   targ = fix(rand*100)+1;
    comI2(iProps(o).PixelIdxList) = cmap(targ,2)*256;
    comI3(iProps(o).PixelIdxList) = cmap(targ,3)*256;    
end
comI(:,:,:,1) = colI(:,:,:,1);
comI(:,:,:,2) = comI2;
comI(:,:,:,3) = comI3;



%%

% 
% for i = 1:size(labI,3)
%    image(uint8(squeeze(comI(:,:,i,:)))),pause(.1) 
% end
% %



for i = 1:size(labI,3)
  imwrite(uint8(squeeze(comI(:,:,i,:))),[TPN iNam{i}]),pause(.1) 
end
