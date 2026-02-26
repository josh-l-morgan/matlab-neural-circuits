if 1
    
startDir = 'E:\LGNs1\exampleRawImages\ExampleCollection02\'
[SFN SPN ] = uigetfile([startDir '*.*'])

end

%%
doMed = 0;
doFilt = 1;
doScale = 0;

medKern = [2 1];
fKern = ones(2,1);
fKern = fKern ./ sum(fKern(:));

rescale = [6/4 1];


fileType = 'tif';%'jpg';%'jpg2000';
compression= 5;
quality = 50;

 
winSize = 200;
%% read
colormap gray(256)
I = imread([SPN SFN]);
[ys xs] = size(I);
If = I;


%% filter


if doMed
    If = medfilt2(If,medKern);
end
if doFilt
    If = imfilter(If,fKern);
end
if doScale
    
   If = imresize(If,[round(ys*rescale(1)) round(xs*rescale(2))],'bicubic'); 
    
end


%% write
if strcmp(fileType,'tif');
    TFN = [SFN(1:end-4) '_Filt01.tif'];
    imwrite(If,[SPN TFN])
elseif strcmp(fileType,'jpg')
    TFN = [SFN(1:end-4) '_Filt01.jpg'];
    imwrite(If,[SPN TFN],'quality',quality)
elseif strcmp(fileType,'jpg2000')
    TFN = [SFN(1:end-4) '_Filt01.j2k'];
    imwrite(If,[SPN TFN],'CompressionRatio',compression)
end

%% show result
if 0 

for r = 1:10
    for yc = winSize/2:winSize:ys-winSize/2
        for xc = winSize/2:winSize/2:xs-winSize/2
            window = [round(yc-winSize/2) round(yc+winSize/2); round(xc-winSize/2) round(xc+winSize/2)];
            window(window<1) = 1;
            window(1,2) = min(window(1,2),ys);
            window(2,2) = min(window(2,2),xs);
            
            
            subplot(2,1,1)
            image(I(window(1,1):window(1,2),window(2,1):window(2,2)));
            
            subplot(2,1,2)
            image(If(window(1,1):window(1,2),window(2,1):window(2,2)))
            
            pause;
        end
    end
end

else
    
    yc = ys/2 +100;
    xc = xs/2 -100;
    window = [round(yc-winSize/2) round(yc+winSize/2); round(xc-winSize/2) round(xc+winSize/2)];
    window(window<1) = 1;
    window(1,2) = min(window(1,2),ys);
    window(2,2) = min(window(2,2),xs);
    
    
    subplot(2,1,1)
    image(I(window(1,1):window(1,2),window(2,1):window(2,2)));
    
    subplot(2,1,2)
    image(If(window(1,1):window(1,2),window(2,1):window(2,2)))
    
    pause(.01);
    
    
end





