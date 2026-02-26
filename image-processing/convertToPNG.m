SPN = 'F:\gxI_13-15\Aligned Stack\'
TPN = 'F:\gxI_13-15\Aligned Stack CroppedJM\'

dSPN = dir([SPN '*.tif']);
inams = {dSPN(:).name};
crop = [720 11312; 3488 12304];

dify = crop(1,2)-crop(1,1)
crop(1,2) = crop(1,1) + ceil(dify/64)*64 -1;

difx = crop(2,2)-crop(2,1)
crop(2,2) = crop(2,1) + ceil(difx/64)*64 -1;

for i = 1:length(inams)
   nam = inams{i};
    I = imread([SPN nam]);
    I = I(crop(1,1):crop(1,2), crop(2,1):crop(2,2));
    filename = [TPN nam(1:end-3) 'png'];
    imwrite(I,filename);
    
end
