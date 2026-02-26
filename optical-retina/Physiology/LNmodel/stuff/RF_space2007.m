STA = spikes(i).STA;
STA(1,:) = 0.5;
tpeak = find (std(STA) == max(std(STA)));
fpeak = mean(STA(:,tpeak-2:tpeak+2)');


if ( mean(mean(STA))- min(min(STA))  <  max(max(STA))- mean(mean(STA)) )        %FOR ON CELLS
    RF_ON = find ( fpeak >= ( mean(fpeak) + 3*std(fpeak) ) );

    absCent = find (fpeak == max(fpeak));
    mask = zeros(80,60);
    mask(absCent) = 1;
    SE = strel ('square', 8);
    ROI = imdilate(mask,SE);
    VROI = find (ROI > 0);
    T = [-800:10];


    Cent = zeros(80,60);
    Cent(RF_ON) = 1;
    SE1 = strel ('diamond',1);
    Int = imdilate(Cent,SE1);
    SE2 = strel ('diamond',4);
    Sur = imdilate(Cent,SE2);
    Sur = Sur - Int;
    RF_OFF = find (Sur > 0);


    area = (66)^2*length(RF_ON);
    diameter = sqrt(area/pi);

    figure(1)
    subplot (3,1,1)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([min(min(STA)), max(max(STA))]);

    subplot (3,1,2)
    abstrus = .5 + zeros(80,60);
    abstrus (RF_ON) = 1;
    abstrus (RF_OFF) = 0;
    pcolor (abstrus')
    shading flat
    colormap gray
    caxis([0 1]);

    subplot (3,1,3)
    pcolor (ROI')
    shading flat
    colormap gray
    caxis([0 1]);

    figure(2)
    for i=1:length(VROI)
        pix(i,:) = STA(VROI(i),:);
        A = zeros(80,60);
        A(VROI(i)) = 1;
        C = zeros(80,60);
        C (RF_ON) =1;
        S = zeros(80,60);
        S (RF_OFF) =1;
        subplot (8,8,i)
        if any(any(A.*C > 0))
            plot(T,pix(i,:),'r');
        elseif any(any(A.*S > 0))
            plot(T,pix(i,:),'b');
        else
            plot(T,pix(i,:),'k');
        end
        axis([-800 10 min(min(STA)) max(max(STA))]);
%         axis tight
        axis off
    end
    
    %temporal receptive field averages
    for j=1:length(RF_ON)
        ON(j,:) = STA(RF_ON(j),:);
    end
    %
    for j=1:length(RF_OFF)
        OFF(j,:) = STA(RF_OFF(j),:);
    end
    %
    T_ON = mean(ON);
    TTP_C = 800-find(T_ON==(max(T_ON)));
    TTP_C_ant = 800-find(T_ON==(min(T_ON)));
    T_OFF = mean(OFF);
    TTP_S = 800-find(T_OFF==(min(T_OFF)));
    TTP_S_ant = 800-find(T_OFF==(max(T_OFF)));
    T = [-800:10];
    %
    figure(3)
    subplot (3,1,1)
    plot(T,T_ON,'r')
    axis tight
    %
    subplot (3,1,2)
    plot(T,T_OFF,'b')
    axis tight

    subplot (3,1,3)
    plot(T,T_ON,'r')
    hold
    plot(T,T_OFF,'b')
    axis tight



elseif ( mean(mean(STA))- min(min(STA))  >  max(max(STA))- mean(mean(STA)) )      %FOR OFF CELLS
    RF_OFF = find ( fpeak <= ( mean(fpeak) - 3*std(fpeak) ) );

    absCent = find (fpeak == min(fpeak));
    mask = zeros(80,60);
    mask(absCent) = 1;
    SE = strel ('square', 15);
    ROI = imdilate(mask,SE);
    VROI = find (ROI > 0);
    T = [-800:10];


    Cent = zeros(80,60);
    Cent(RF_OFF) = 1;
    SE1 = strel ('diamond',1);
    Int = imdilate(Cent,SE1);
    SE2 = strel ('diamond',4);
    Sur = imdilate(Cent,SE2);
    Sur = Sur - Int;
    RF_ON = find (Sur > 0);
    
    area=(66)^2*length(RF_OFF);
    diameter = sqrt(area/pi);

    figure(1)
    subplot (3,1,1)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([min(min(STA)), max(max(STA))]);

    subplot (3,1,2)
    abstrus = .5 + zeros(80,60);
    abstrus (RF_ON) = 1;
    abstrus (RF_OFF) = 0;
    pcolor (abstrus')
    shading flat
    colormap gray
    caxis([0 1]);

    subplot (3,1,3)
    pcolor ((-ROI+1)')
    shading flat
    colormap gray
    caxis([0 1]);

    figure(2)
    for i=1:length(VROI)
        pix(i,:) = STA(VROI(i),:);
        A = zeros(80,60);
        A(VROI(i)) = 1;
        C = zeros(80,60);
        C (RF_OFF) =1;
        S = zeros(80,60);
        S (RF_ON) =1;
        subplot (15,15,i)
        if any(any(A.*C > 0))
            plot(T,pix(i,:),'r');
        elseif any(any(A.*S > 0))
            plot(T,pix(i,:),'b');
        else
            plot(T,pix(i,:),'k');
        end
        axis([-800 10 min(min(STA)) max(max(STA))]);
%         axis tight        
        axis off
    end
    
    
    %temporal receptive field averages
    for j=1:length(RF_ON)
        ON(j,:) = STA(RF_ON(j),:);
    end
    %
    for j=1:length(RF_OFF)
        OFF(j,:) = STA(RF_OFF(j),:);
    end
    %
    T_ON = mean(ON);
    TTP_S = 800-find(T_ON==(max(T_ON)));
    TTP_S_ant = 800-find(T_ON==(min(T_ON)));
    T_OFF = mean(OFF);
    TTP_C = 800-find(T_OFF==(min(T_OFF)));
    TTP_C_ant = 800-find(T_OFF==(max(T_OFF)));
    T = [-800:10];
    %
    figure(3)
    subplot (3,1,1)
    plot(T,T_OFF,'r')
    axis tight
    %
    subplot (3,1,2)
    plot(T,T_ON,'b')
    axis tight

    subplot (3,1,3)
    plot(T,T_OFF,'r')
    hold
    plot(T,T_ON,'b')
    axis tight

else
    disp(error)
end


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