function[] = tt_filtStack(app,filtType,kSize)
global globTT

Ic = globTT.I.tab{globTT.active.ID};
chan = find([app.twkR.Value app.twkG.Value app.twkB.Value]);

%Switch to filter tab
id = globTT.active.ID;
id = globTT.lu.filt(id);
app.mainTab.SelectedTab = app.mainTab.Children(id);
pickTabP(app,id)

%% run filter
colormap(globTT.active.ax,gray(256));
for c = chan
    
    I = squeeze(Ic(:,:,c,:));
    
    if strcmp(filtType,'median')
        if kSize(3) <2; %2D filter
            for p = 1:size(I,3);
                I(:,:,p) = medfilt2(I(:,:,p),[kSize(1) kSize(2)]);
                image(globTT.active.ax,uint8(I(:,:,p)));
                pause(.001)
            end
        else
            I = medfilt3(I,kSize);
        end
        
    elseif strcmp(filtType,'gaussian')
        if kSize(3) 
            I = imgaussfilt3(I,kSize(1));
        else
            for p = 1:size(I,3);
                I(:,:,p) = imgaussfilt(I(:,:,p),kSize(1));
                 image(globTT.active.ax,uint8(I(:,:,p)));
                pause(.001)
            end
            
        end
        
    elseif strcmp(filtType,'mexhat')
      oldMax = max(I(:));
        if kSize(3)
            I1 = imgaussfilt3(I,kSize(1));
            I2 = imgaussfilt3(I,kSize(2));
            Id = I1 - I2;
            I = Id;
        else
            for p = 1:size(I,3)
                I1 = imgaussfilt(I(:,:,p),kSize(1));
                I2 = imgaussfilt(I(:,:,p),kSize(2));
                I(:,:,p) = I1-I2;
                maxP = max(max(I(:,:,p)));
                image(globTT.active.ax,uint8(I(:,:,p)*256/maxP));
                pause(.001)
            end
        end
       % I = I - min(I(:));
        I(I<0) = 0;
        newMax = max(I(:));
        I = I  * 256/newMax;
        
    end
    
    Ic(:,:,c,:) = I;
    
end %chan


globTT.I.tab{id} = Ic;






