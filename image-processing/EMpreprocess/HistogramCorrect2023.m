%%Correct contrast for vsvi mipmaps



%% Define folders
SPN = 'Z:\Active\morganLab\DATA\jxQ_PMC_forPhil\'

mipSourceFold = 'mipmaps';
mipTargetFold = 'mipmapsBC';
contrastMipLevel = '4';
imageType = 'png';

TPN = [SPN mipTargetFold];
if ~exist(TPN),mkdir(TPN);end

%% Set variables
invert = 1;
topVal = 210;
botVal = 80;
topRat = .98;
botRat = .02;

%% Gather data
secDir = dir([SPN mipSourceFold]);
secNames = {secDir(3:end).name};

hRange = 0:255;
bVals = zeros(length(secNames),1);
tVals = bVals;
for i = 1:length(secNames)    
    disp(sprintf('%d of %d',i,length(secNames)))
    trackVals = hRange * 0;
    tDir = [SPN mipSourceFold '\' secNames{i} '\' contrastMipLevel '\'];
    iDir = dir([tDir '*.' imageType]);
    iNams = {iDir.name};
    for n = 1:length(iNams)
        I = imread([tDir iNams{n}]);
        trackVals = trackVals + hist(I(:),hRange); 
    end
    trackVals(1) = 0;
    trackVals(end) = 0;
    cSum = cumsum(trackVals);
    cSum = cSum/max(cSum);

    bSum = find(cSum<=botRat);
    bVals(i) = bSum(end);
    tSum = find(cSum>=topRat);
    tVals(i) = tSum(1);

end
colormap gray(256)

%% Adjust values

for i = 1:length(secNames)
    

    topNow = tVals(i);
    botNow = bVals(i);
    
    Ifd2 = Ifd2 - botNow;
    Ifd2 = Ifd2 * (topVal - botVal)/(topNow - botNow);
    Ifd2 = Ifd2 + botVal;
    
    %%image
    image(Ifd2);
    pause(.01)
end



%% Old stuff
inam = {iDir.name};
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