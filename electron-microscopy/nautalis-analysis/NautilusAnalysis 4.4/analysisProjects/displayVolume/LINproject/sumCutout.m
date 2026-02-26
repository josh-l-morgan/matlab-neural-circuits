

%%

SPN = 'D:\LGNs1\Segmentation\VAST\S4\joshm\export\export_16+04+07_320\'


idir = dir([SPN '*.png']);

inams = {idir.name}

clear I
for i = 1:length(inams)
    
   I(:,:,:,i) = double(imread([SPN inams{i}])); 
    
end

Isum = sum(I,4);
image(uint8(Isum/40))






