

colormap gray(256)

SPN = 'E:\Pratyush\Processing\testProcessW12\'
TPN = 'E:\Pratyush\Processing\testProcessW12_medDown\'
if ~exist(TPN),mkdir(TPN);end

iDir = dir([SPN '*.tif'])
inam = {iDir.name};
invert = 1;
downSampFactor = 4;

iInfo = imfinfo([SPN inam{1}]);
ys = iInfo.Height;
xs = iInfo.Width;

secNum = length(inam);

for i = 1:secNum
    sprintf('%d of %d',i,secNum)
     filename = sprintf('%s%05.0f_ds%d.png',TPN,i,downSampFactor);
    if 1;%~exist(filename,'file')
    I = imread([SPN inam{i}]);
    
    %% 
    If = single(I);
    %If = If(12200:12600,15200:15600);
    
    if invert
        If = 255-If;
    end
    
    %%median filter
    if 0
        If = medfilt2(If,[3 3]);
    else
        gkern = fspecial('gaussian',10,1);
        If = imfilter(If,gkern,'same');
    end
    
    Ifd = imresize(If,1/4,'bicubic');
    image(Ifd)
    
    
    bkern = fspecial('gaussian',1000,250);
    IfB = imfilter(Ifd,bkern,'same','symmetric');
    Ifd2 = Ifd-IfB;
    image(Ifd2+100)
    
    %%histogram correction
    topVal = 210;
    botVal = 80;
    topRat = .02;
    botRat = .98;
    val = sort(Ifd2(:),'descend');
    L = length(val);
    topNow = val(round(topRat * L));
    botNow = val(round(botRat * L));
    
    Ifd2 = Ifd2 - botNow;
    Ifd2 = Ifd2 * (topVal - botVal)/(topNow - botNow);
    Ifd2 = Ifd2 + botVal;
    
    %%image
    image(Ifd2);
    pause(.01)
        
        
    %% save
    imwrite(uint8(Ifd2),filename)
    end
end