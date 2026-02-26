SPN = 'D:\Human\HuB\JennaFil_HuB5\ProbableRGC.oif.files\'
TPN = [SPN(1:end-1) '_stack\']
if ~exist(TPN),mkdir(TPN);end
sDir = dir([SPN '*.tif']);

nams = {sDir.name};


for i = 1:length(nams)
   nam = nams{i}
    c = regexp(nam,'C');
    z = regexp(nam,'Z');
    dotT = regexp(nam,'.tif');
    col(i) = str2num(nam(c+1:z-1));
    plane(i) = str2num(nam(z+1:dotT-1));
end

for i = 1:max(plane)
    
    R = find((col == 1) & (plane == i));
    if ~isempty(R)
        I = imread([SPN nams{R}]);
        if strcmp(class(I),'uint16')
            Ic(:,:,1) = I/256;
        else
            Ic(:,:,1) = I;
        end
    end
    
     G = find((col == 2) & (plane == i));
    if ~isempty(G)
        I = imread([SPN nams{G}]);
        if strcmp(class(I),'uint16')
            Ic(:,:,2) = I/16;
        else
            Ic(:,:,2) = I;
        end
    end
    
     B = find((col == 3) & (plane == i));
    if ~isempty(B)
        I = imread([SPN nams{B}]);
        if strcmp(class(I),'uint16')
            Ic(:,:,3) = I/16;
        else
            Ic(:,:,3) = I;
        end
    end
    
    Ic = uint8(Ic);
    %image(Ic*10),pause(.01)
    newName = sprintf('plane_%04.0f.tif',i);
    imwrite(Ic,[TPN newName]);
    
end





    