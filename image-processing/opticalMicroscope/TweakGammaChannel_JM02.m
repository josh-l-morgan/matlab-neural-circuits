clear all
SPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack2\'
TPN = 'E:\Pratyush\Optical\confocal\th_tdt_121019_right_trap\20X_slow\03_oifStack_tweak6\'
mkdir(TPN) 

dSPN = dir([SPN '*.tif']);
inams = {dSPN.name};
p = length(inams);

shiftHisto = 0; %Use rolling histogram correction (1) or single correction for whole stack

gamma1 = .8;
gamma2 = .8;
topVal = 220; %(assuming 8 bit)
botVal = 20;
topRat = .001;
botRat = .999;

topBit = 2^16 - 1;
topVal = topVal/255*topBit


binRad = 30;

info = imfinfo([SPN inams{1}]);

I1 = zeros(info.Height,info.Width,p);
I2 = I1;


for i = 1:p
    I = double(imread([SPN inams{i}]));
    I1(:,:,i) = medfilt2(I(:,:,1),[3 3]);
    I2(:,:,i) = medfilt2(I(:,:,2),[3 3]);

end

stackMin1 = min(I1(:));
stackMin2 = min(I2(:));

stackMax1 = max(I1(:));
stackMax2 = max(I2(:));


for i = 1:p
    disp(sprintf('%d of %d',i,p))
    start = max(1,i-binRad);
    stop = min(p,i+binRad);
    
    %%%%%%%%%%%%%% Channel 1
    
    Is = I1(:,:,start:stop);
    Ic = I1(:,:,i);
    
    %%histogram correction
    if shiftHisto
        val = sort(Is(:),'descend');
        L = length(val);
        topNow = val(round(topRat * L));
        botNow = val(round(botRat * L));
    else
        topNow = stackMax1;
        botNow = stackMin1;
    end
    
    Ic = Ic - botNow;
    Ic = Ic * (topVal - botVal)/(topNow - botNow);
    Ic = Ic + botVal;
    
    image(uint8(Ic))
    
    Ig = Ic;
    Ig(Ig<1) = 1;
    Ig(Ig>topBit) = topBit;
    Ig = (Ig/topBit).^gamma1*topBit;
    image(uint8(Ig))
    
    Icomb(:,:,1) = Ig;
    
    %%%%%%%%% Channel 2
    Is = I2(:,:,start:stop);
    Ic = I2(:,:,i);
    
    %%histogram correction
    if shiftHisto
        val = sort(Is(:),'descend');
        L = length(val);
        topNow = val(round(topRat * L));
        botNow = val(round(botRat * L));
    else
        topNow = stackMax2;
        botNow = stackMin2;
    end
    
    Ic = Ic - botNow;
    Ic = Ic * (topVal - botVal)/(topNow - botNow);
    Ic = Ic + botVal;
    
    image(uint8(Ic))
    
    Ig = Ic;
    Ig(Ig<1) = 1;
    Ig(Ig>topBit) = topBit;
    Ig = (Ig/topBit).^gamma1*topBit;
    image(uint8(Ig))
    
    Icomb(:,:,2) = Ig;
    
    Icomb(:,:,3) = Ig * 3;
    pause(.01)
    filename = sprintf('%s%05.0f.tif',TPN,i);
    imwrite(uint16(Icomb),filename);
    
end





