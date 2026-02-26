clear all
lastDir = 'E:\SeanMcCracken\50um_HighRes\'
[SFN SPN] = uigetfile([lastDir '*']);

r = bfopen([SPN SFN]);


%%Parse labels
imageCell = r{1};
clear secIDs colIDs
for i = 1:size(imageCell,1)
    nam = imageCell{i,2};
    
    C = regexp(nam,'C=');
    colIDs(i) = str2num(nam(C(end)+2));
    
    Z = regexp(nam,'Z=');
    Zstring = nam(Z(end)+2:end);
    Zslash = regexp(Zstring,'/');
    secIDs(i) = str2num(Zstring(1:Zslash(1)-1));    
end
secs = unique(secIDs);
cols = unique(colIDs);

%%Make  stack

if length(cols)>3
    useCol = [1 3 4];
elseif length(cols) ==3
    useCol = cols;
elseif length(cols) == 2
    useCol = [cols(1) cols(2) cols(1)];
    
else
    useCol = [cols cols cols];
end


I1 = imageCell{1,1};
clear I
if 1
    Istack = zeros(size(I1,1),size(I1,2),3,length(secs));
    
    for i = 1:length(secs)
        isSec = find(secIDs==secs(i));
        secCol = colIDs(isSec);
        
        for c = 1: length(useCol)
            
            targCol = find(secCol==useCol(c),1);
            p = isSec(targCol);
            I(:,:,c) = imageCell{p,1};
            
        end
        Istack(:,:,:,i) = I;
        
        image(I),drawnow
               
    end
end

%Istack = Istack(:,:,:,1:35);

if max(Istack(:))>255
    scaleI = 1/2^8;
else
    scaleI = 1;
end





%%Med filter
filtChan = [1 2 3];
Ifilt = Istack;
for i = 1:size(Ifilt,4);
    i
    for c = 1:length(filtChan)
        Ifilt(:,:,filtChan(c),i) = medfilt2(Ifilt(:,:,filtChan(c),i),[2 2]);
    end  
end

Imax = max(Ifilt,[],4) ;
for c = 1:3
    Imax(:,:,c) = Imax(:,:,c) * 255/max(max(Imax(:,:,c)));
end
image(uint8(Imax)), drawnow
maxName = ['Max_' SFN(1:end-4) '.tif'];
imwrite(uint8(Imax),[SPN maxName]);



%%Write
Iwrite  = Ifilt;
TPN = [SPN SFN(1:end-4) '_stack\'];
mkdir(TPN)
numZ = size(Istack,4);
for i = 1:numZ;
    i
            fileName = sprintf('%03.0f.tif',i);
            
            Is = Iwrite(:,:,:,i) * scaleI;
           imwrite(uint8(Is),[TPN fileName]);
end

doLoop = 1;

if doLoop
    
    for i = 1:numZ;
        i
        fileName = sprintf('%03.0f.tif',numZ+i);
        
        Is = Iwrite(:,:,:,numZ-i+1) * scaleI;
        imwrite(uint8(Is),[TPN fileName]);
    end
    
end









