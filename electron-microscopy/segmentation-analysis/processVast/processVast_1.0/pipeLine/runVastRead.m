
SPN = 'D:\LGNs1\segmentation\VAST\S8\joshm\export_InsideOut\';
TPN = 'D:\LGNs1\segmentation\VAST\S8\joshm\export_InsideOutMat\';
if ~exist(TPN,'dir'),mkdir(TPN),end


%% creat obj struct from tiles
if exist([TPN 'obj.mat']);
    
    load([TPN 'obj.mat']);

else
    obj.subs = vastTiles2subs(SPN);
    
    getText = dir([SPN '*.txt']);
    if ~isempty(getText)
        obj.ids = readVastColors([SPN getText(1).name]);
    end
    obj.images = viewObj01(obj,1);  
    save([TPN 'obj.mat'],'obj','-v7.3')
    
end

image(uint8(obj.images.colorOb))
imwrite(uint8(obj.images.colorOb),[TPN 'colorSegments.tif']); 

