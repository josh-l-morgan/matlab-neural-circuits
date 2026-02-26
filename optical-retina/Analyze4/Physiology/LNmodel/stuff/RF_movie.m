% receptive field movie
clear STA h* low mg M MovF 
STA = spikes(i).STA;
STA(1,:)=0.5;
w = 811;
mg = mean(STA(:));
low = min(STA(:));
high = max(STA(:));
hFig1 = figure('Position',[50 100 500 500]);
if (mg-low) >= (high-mg)
    for j=1:w
        MovF = (reshape(STA(:,j),80,60))';
        pcolor(MovF)
        shading flat
        colormap gray
        caxis([low, 1 - low]);
%         axis([26 54 20 48])
        axis off
        M(j) = getframe;
    end
else
    for j=1:w
        MovF = (reshape(STA(:,j),80,60))';
        pcolor(MovF)
        shading flat
        colormap gray
        caxis([1 - high, high]);
%         axis([26 54 20 48])
        axis off
        M(j) = getframe;
    end
end
