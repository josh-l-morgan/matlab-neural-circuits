%this creates plot of the temporal structure of the receptive field center
%first-pixel anomaly correction
STA = spikes(i).STA;
timeWindow = size(STA);
meanGray = 0.5;
STA(1,:) = meanGray*ones(1,timeWindow(2));

tpeak = find (std(STA) == max(std(STA)));
fpeak = mean(STA(:,tpeak-2:tpeak+2)');

mg = mean(STA(:));
low = min(STA(:));
high = max(STA(:));


if (high-mg) > (mg-low)        %FOR ON CELLS
    RF_ON = find ( fpeak >= ( mean(fpeak) + 4*std(fpeak) ) );

    %spatial structure at temporal max
    figure(1)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([1 - high, high]);
%     axis([26 46 20 40])
    axis off




    %temporal receptive field averages
    for j=1:length(RF_ON)
        ON(j,:) = STA(RF_ON(j),:);
    end

    T_ON = mean(ON)-.5;
    TTP_C = 800-find(T_ON==(max(T_ON)));
    TTP_C_ant = 800-find(T_ON==(min(T_ON)));
    T = [-800:10];

    figure(2)
    plot(T,T_ON/max(T_ON),'k')
    axis tight
    %



elseif (mg-low) >= (high-mg)      %FOR OFF CELLS
    RF_OFF = find ( fpeak <= ( mean(fpeak) - 4*std(fpeak) ) );

    %spatial structure at temporal max
    figure(1)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([low, 1 - low]);
%     axis([26 46 20 40])
    axis off



    %temporal receptive field averages
    for j=1:length(RF_OFF)
        OFF(j,:) = STA(RF_OFF(j),:);
    end

    T_OFF = mean(OFF)-.5;
    TTP_C = 800-find(T_OFF==(min(T_OFF)));
    TTP_C_ant = 800-find(T_OFF==(max(T_OFF)));
    T = [-800:10];
    %
    figure(2)
    plot(T,T_OFF/abs(min(T_OFF)),'k')
    axis tight

else
    disp(error)
end

disp(spikes(i).channel)

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