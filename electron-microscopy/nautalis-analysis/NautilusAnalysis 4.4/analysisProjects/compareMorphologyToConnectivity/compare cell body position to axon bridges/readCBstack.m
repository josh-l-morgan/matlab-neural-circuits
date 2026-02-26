
%SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\exportAutoCB_3editNewCol2\'
%SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\export_CBedit-50um\'
SPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\autoCB_3edit3f_allChecked_export\'
TPN = 'D:\LGNs1\Segmentation\VAST\S8\S8_ds16\CellBodies\matOut_AutoCB\';

obI = makeOBI(SPN,TPN)


onam = obI.nameProps.names;
for i = 1:length(onam)
    isSeg(i) = sum(regexp(onam{i},'Segment'));
end


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

%[labVol numLab] = bwlabeln(idVol>0);
cbObj = bwconncomp(idVol);
props = regionprops(cbObj)
Areas = [props.Area];
FailSize = find(Areas<50);
PassSize = find(Areas>50);

Pass = PassSize & isSeg;
filtVol = idVol*0;
for i = 1:length(Pass)
       filtVol(cbObj.PixelIdxList{Pass(i)}) = i;   
end




%[labVol numLab] = bwlabeln(filtVol>0);

cbObj = bwconncomp(filtVol);
cbObj.vol = filtVol;

save([TPN 'cbObj.mat'],'cbObj','-v7.3')
% 
% for i = 1:cbObj.NumObjects
%     i
%    vals = idVol(cbObj.PixelIdxList{i});
%    isDif(i) = sum(abs(vals(1:end-1)-vals(2:end)));
%     
%    
% end







