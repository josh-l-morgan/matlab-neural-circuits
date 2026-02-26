% receptive field movie
clear STA h* low mg M MovF col row
clc

STA = spikes(i).STA;
STA(1,:)=0.5;
w = 811;
mg = mean(STA(:));
low = min(STA(:));
high = max(STA(:));
tpeak = find (std(STA) == max(std(STA)));
fpeak = mean(STA(:,tpeak-2:tpeak+2),2);
frame = reshape(fpeak,80,60)';

hFig1 = figure('Position',[400 500 500 500]);
if (mg-low) >= (high-mg) %OFF cells
    [row col] = find(frame == min(frame(:)));
    nPixels = length(find(fpeak <= (mean(STA(:)) - 3*std(fpeak(:)))));
    for j=1:500
        MovF = (reshape(STA(:,j+299),80,60))';
        pcolor(MovF)
        shading flat
        colormap gray
        caxis([low, 1 - low]);
        axis([col-8 col+8 row-8 row+8])
        axis off
        M(j) = getframe;
    end
else %ON cells
    [row col] = find(frame == max(frame(:)));
    nPixels = length(find(fpeak >= (mean(STA(:)) + 3*std(fpeak(:)))));
    for j=1:500
        MovF = (reshape(STA(:,j+299),80,60))';
        pcolor(MovF)
        shading flat
        colormap gray
        caxis([1 - high, high]);
        axis([col-8 col+8 row-8 row+8])
        axis off
        M(j) = getframe;
    end
end
% movie(M)
% movie2avi(M,'testbild','compression','none','fps',30)
disp(spikes(i).channel)
% disp(nPixels)
