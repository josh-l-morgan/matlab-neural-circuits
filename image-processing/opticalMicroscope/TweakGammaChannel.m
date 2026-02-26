
SPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack\'
TPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack_tweak2\'
mkdir(TPN)

dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};
p = length(inams);

gamma1 = .8;
gamma2 = .8;

binRad = 1

info = imfinfo([SPN inams{1}]);

I1 = zeros(info.Height,info.Width,p);
I2 = I1;

for i = 1:p
    I = imread([SPN inams{i}]);
    I1(:,:,i) = medfilt2(I(:,:,1),[5 5]);
    I2(:,:,i) = medfilt2(I(:,:,2),[5 5]);
end
    

for i = 1:p
    disp(sprintf('%d of %d',i,p))
    start = max(1,i-binRad);
    stop = min(p,i+binRad);
    
    Is = I1(:,:,start:stop);
    Ic = I1(:,:,i);
    
    %%histogram correction
    topVal = 250;
    botVal = 10;
    topRat = .01;
    botRat = .98;
    val = sort(Is(:),'descend');
    L = length(val);
    topNow = val(round(topRat * L));
    botNow = val(round(botRat * L));
    
    Ic = Ic - botNow;
    Ic = Ic * (topVal - botVal)/(topNow - botNow);
    Ic = Ic + botVal;
    
    image(uint8(Ic))
    
    Ig = Ic;
    Ig(Ig<1) = 1;
    Ig(Ig>255) = 255;
    Ig = (Ig/255).^gamma1*256;
    image(uint8(Ig))
    
    Icomb(:,:,1) = Ig;
    
    %%%%%%%%% Channel 2
     Is = I2(:,:,start:stop);
    Ic = I2(:,:,i);
    
    %%histogram correction
    topVal = 250;
    botVal = 10;
    topRat = .01;
    botRat = .98;
    val = sort(Is(:),'descend');
    L = length(val);
    topNow = val(round(topRat * L));
    botNow = val(round(botRat * L));
    
    Ic = Ic - botNow;
    Ic = Ic * (topVal - botVal)/(topNow - botNow);
    Ic = Ic + botVal;
    
    image(uint8(Ic))
    
    Ig = Ic;
    Ig(Ig<1) = 1;
    Ig(Ig>255) = 255;
    Ig = (Ig/255).^gamma1*256;
    image(uint8(Ig))
    
    Icomb(:,:,2) = Ig;
    
    Icomb(:,:,3) = Ig * 5;
 pause(.01)
    filename = sprintf('%s%05.0f.tif',TPN,i);
    imwrite(uint8(Icomb),filename);
    
end





