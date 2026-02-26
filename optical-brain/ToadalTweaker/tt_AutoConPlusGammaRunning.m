function[] = tt_AutoConPlusGammaRunning(app)

global globTT

shouldRam = 1;
shouldWrite = 0;
colormap(globTT.active.ax,gray(256));



I = globTT.I.tab{globTT.active.ID};
[ys xs cs p] = size(I);
if globTT.active.ID<4
    targID = 2;
else
    targID = 5;
end

app.mainTab.SelectedTab = app.mainTab.Children(targID);
globTT.active.ID = targID;
globTT.active.ax = app.mainTab.SelectedTab.Children(1);

binRad = app.edit_autoConBin.Value;
shiftHisto = binRad>0;%Use rolling histogram correction (1) or single correction for whole stack

gamma1 = 1;
topVal = 220; %(assuming 8 bit)
botVal = 20;
topRat = .001;
botRat = .999;
% 
% topBit = 2^16 - 1;
% topVal = topVal/255*topBit

biggest = max(I(:));
for i = 3:20
    bitCeiling = 2^i -1;
    if bitCeiling >=biggest 
        break
    end
end


%info = imfinfo([SPN inams{1}]);
Ia = I;
for c = 1:cs
 
    I1 = squeeze(I(:,:,c,:));
    stackMax = max(I1(:));
    stackMin = min(I1(:));
    
    for i = 1:p
        disp(sprintf('%d of %d',i,p))
        start = max(1,i-binRad);
        stop = min(p,i+binRad);
        
        %%%%%%%%%%%%%% Channel 1
        
        Is = squeeze(I1(:,:,start:stop));
        Ic = I1(:,:,i);
        
        %%histogram correction
        if shiftHisto
            val = sort(Is(:),'descend');
            L = length(val);
            topNow = val(round(topRat * L));
            botNow = val(round(botRat * L));
        else
            topNow = stackMax;
            botNow = stackMin;
        end
        
        Ic = Ic - botNow;
        Ic = Ic * (topVal - botVal)/(topNow - botNow);
        Ic = Ic + botVal;
        
        
        Ig = Ic;
        Ig(Ig<0) = 0;
        %Ig(Ig>topBit) = topBit;
        oldTop = max(Ig(:));
        Ig = Ig.^gamma1;
        Ig = Ig * oldTop/max(Ig(:));
       % image(globTT.active.ax,uint8(Ig*256/topVal))
        image(globTT.active.ax,uint8(Ig))
        Ia(:,:,c,i) = Ig;
        
        pause(.01)
        if shouldWrite
            filename = sprintf('%s%05.0f.tif',TPN,i);
            imwrite(uint16(Ia),filename);
        end
        
    end
    
end

globTT.I.tab{targID} = Ia;
    
    
