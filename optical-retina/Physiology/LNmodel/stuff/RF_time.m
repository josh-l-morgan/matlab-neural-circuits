% function [timeKernel] = receptivefieldtime(STA)


peakTime = find (std(STA) == max(std(STA)));
peakFrame = mean(STA(:,peakTime-1:peakTime+1),2);


if ( mean(mean(STA))- min(min(STA))  <  max(max(STA))- mean(mean(STA)) )        %FOR ON CELLS
    RF_ON = find ( fpeak >= ( mean(fpeak) + 3*std(fpeak) ) );

    %spatial structure at temporal max
    figure(2)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([min(min(STA)), max(max(STA))]);


    
    %temporal receptive field averages
    for j=1:length(RF_ON)
        ON(j,:) = STA(RF_ON(j),:);
    end
    
    T_ON = mean(ON)-.5;
    TTP_C = 800-find(T_ON==(max(T_ON)));
    TTP_C_ant = 800-find(T_ON==(min(T_ON)));
    T = [-800:10];
    
    figure(1)
    plot(T,T_ON,'m')
    hold on
    axis tight
    %



elseif ( mean(mean(STA))- min(min(STA))  >  max(max(STA))- mean(mean(STA)) )      %FOR OFF CELLS
    RF_OFF = find ( fpeak <= ( mean(fpeak) - 3*std(fpeak) ) );

    %spatial structure at temporal max
    figure(2)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([min(min(STA)), max(max(STA))]);

    
    
    %temporal receptive field averages
    for j=1:length(RF_OFF)
        OFF(j,:) = STA(RF_OFF(j),:);
    end
    
    T_OFF = mean(OFF)-.5;
    TTP_C = 800-find(T_OFF==(min(T_OFF)));
    TTP_C_ant = 800-find(T_OFF==(max(T_OFF)));
    T = [-800:10];
    %
    figure(1)
    plot(T,T_OFF,'m')
    hold on
    axis tight

else
    disp(error)
end

clear T* O* RF* snap
% %RF map
%
% % %spatial receptive field
% for j=1:w
%     MovM = (reshape(STA(:,j),80,60))';
%     pcolor(MovM)
%     shading flat
%     colormap gray
%     caxis([min(min(STA)), max(max(STA))]);
%     MovS(j) = getframe;
% end
% % snap = reshape(fpeak, 80, 60)';
% % mesh (snap)
% % shading flat
% % colormap gray
% % caxis([min(min(fpeak)), max(max(fpeak))]);