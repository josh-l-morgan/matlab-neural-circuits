
SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\exportAutoCB_3editNewCol2\'
TPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\matOut_AutoCB\';


imageNames = dir([SPN '*.png']);

testI = imread([SPN imageNames(1).name]);
[ys xs cs] = size(testI);

idVol = zeros(ys, xs, length(imageNames));

for i = 1:length(imageNames)
    i
    Iraw = double(imread([SPN imageNames(i).name]));
    I = Iraw(:,:,1) * 256^2 + Iraw(:,:,2)*256 + Iraw(:,:,3);
    %image(uint8(Iraw)),pause(.01)
    idVol(:,:,i) = I;
end

uLab = unique(idVol);

[labVol numLab] = bwlabeln(idVol>0);
cbObj = bwconncomp(idVol);
cbObj.vol = idVol;

save([TPN 'cbObj.mat'],'cbObj','-v7.3')
% 
% for i = 1:cbObj.NumObjects
%     i
%    vals = idVol(cbObj.PixelIdxList{i});
%    isDif(i) = sum(abs(vals(1:end-1)-vals(2:end)));
%     
%    
% end







