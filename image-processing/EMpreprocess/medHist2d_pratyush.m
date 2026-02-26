

colormap gray(256)

SPN = 'E:\ixD-LGN1\MasterRaw\Waf002_r1-c1\'
TPN = 'E:\ixD-LGN1\MasterRaw\Waf002_r1-c1_med\'
if ~exist(TPN),mkdir(TPN);end

inam = dir([SPN '*.tif'])
invert = 1;

I = imread([SPN inam(1).name]);
[ys xs a] = size(I);

secNum = length(inam);

for i = 1:secNum
    sprintf('%d of %d',i,secNum)
     filename = sprintf('%s%05.0f.png',TPN,i);
    if ~exist(filename,'file')
    I = imread([SPN inam(i).name]);
    
    %% 
    If = I;
    If = If(2200:2600,5200:5600);
    
    if invert
        If = double(255-If);
    end
        
    %%median filter
    If = medfilt2(If,[2 2]);
    
    %%histogram correction
    topVal = 210;
    botVal = 80;
    topRat = .1;
    botRat = .9;
    val = sort(If(:),'descend');
    L = length(val);
    topNow = val(round(topRat * L));
    botNow = val(round(botRat * L));
    
    If = If - botNow;
    If = If * (topVal - botVal)/(topNow - botNow);
    If = If + botVal;
    
    %%image
    image(If);
    pause(.01)
    
        
    %% save
    imwrite(uint8(If),filename)
    end
end