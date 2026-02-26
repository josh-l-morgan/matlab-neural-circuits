CPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\02.oif.files\'
TPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\02_oifStack2\'
mkdir(TPN)

dCPN = dir([CPN '*.tif'])
inams = {dCPN.name}
clear color slice
for i  = 1:length(inams)
   
    nam = inams{i};
    
    c = regexp(nam,'C');
    d = regexp(nam,'.tif');
    z = regexp(nam,'Z');
    
    color(i) = str2num(nam(c+1:z-1));
    slice(i) = str2num(nam(z+1:d-1));
    
    
end

slices = unique(slice);

for i = 1:length(slices)
   
    clear I;
    targ = find((slice == slices(i)) & ( color == 1));
    if ~isempty(targ)
        I = imread([CPN inams{targ}]);
        Ic(:,:,1) = I;
    end
    
    targ = find((slice == slices(i)) & ( color == 2));
    if ~isempty(targ)
        I = imread([CPN inams{targ}]);
        Ic(:,:,2) = I;
    end
    
    targ = find((slice == slices(i)) & ( color == 3));
    if ~isempty(targ)
        I = imread([CPN inams{targ}]);
        Ic(:,:,3) = I;
    end
    if size(Ic,3)<3
        Ic(1,1,3) = 0;
    end
    
    imwrite(uint16(Ic/16),sprintf('%s%05.0f.tif',TPN,slices(i)))
    
end



