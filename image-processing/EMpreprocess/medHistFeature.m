

colormap gray(256)

SPN = 'D:\KarlsRetina\HxR\'
TPN = 'D:\LGN\p7_mouseLGN_med_test\'
showTest = 1;
compressionRatio = 10;
medFiltSize = 2;

if ~exist(TPN),mkdir(TPN);end

inam = dir([SPN '*.tif'])
invert = 1;

tic
'reading raw'

I = imread([SPN inam(1).name]);
toc

[ys xs a] = size(I);

secNum = length(inam);

for i = 1:1%%secNum
    sprintf('%d of %d',i,secNum)
     filename = sprintf('%s%05.0f.j2k',TPN,i);
   
     if 1%~exist(filename,'file')
    I = imread([SPN inam(i).name]);
    
    
    
    %% 
    If = I;
    If = If(10000:12450,10000:12450);
    
    if invert
        If = double(255-If);
    end
    if showTest
        subplot(2,2,1)
        image(If),    pause(.01)
    end    
    %%median filter
    tic
    'med filt'
    If = medfilt2(If,[medFiltSize medFiltSize]);
    toc
    if showTest
        subplot(2,2,2)
        image(If),    pause(.01)
    end    
    
    %%histogram correction
    tic
    'hist correction'
    topVal = 210;
    botVal = 100;
    topRat = .1;
    botRat = .9;
    val = sort(If(:),'descend');
    L = length(val);
    topNow = val(round(topRat * L));
    botNow = val(round(botRat * L));
    
    If = If - botNow;
    If = If * (topVal - botVal)/(topNow - botNow);
    If = If + botVal;
    toc
    if showTest
        subplot(2,2,3)
        image(If),    pause(.01)
    end  
        
        
    %%save
    %imwrite(uint8(If),filename)
    tic
    'writing jp2k'
    imwrite(uint8(If),filename,'CompressionRatio',compressionRatio);
    toc
    
    if showTest
        subplot(2,2,4)
        tic
        'reading jp2k'
        Ij = imread(filename);
        toc
        image(Ij),    pause(.01)
    end
    
    
    Is = 
    
    subplot(2,2,1)
    image(Ij),pause(.01)
    
    regions = detectMSERFeatures(Ij);
    
    subplot(2,2,2)
    image(Ij),hold on;
    plot(regions,'showPixelList',true,'showEllipses',true);pause(.01)

    
    
    
     end
end