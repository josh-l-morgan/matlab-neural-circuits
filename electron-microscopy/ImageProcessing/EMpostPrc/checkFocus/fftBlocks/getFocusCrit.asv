%% Get Image FFT parameters
%{
1) load image
2) grid FFT
3) display FFTs
4) select high and low frequency region for FFT
5) select threshold based on distribution of ratios
%}
clear all
colormap gray(256)

[TFN TPN] = GetMyFile;



sS = 100;  % Define pixel length of area used for FFT sampling
minSigArea = .1; % proportion of image that must have signal



%% Start Reading


    I = double(imread([TPN TFN]));
    %I = I(1:3000,:);
    [ys xs] = size(I);

    iFFTs = [];

    yStep = 1:sS:ys-sS;
    xStep = 1:sS:xs-sS;
    fftMos = zeros(length(yStep),length(xStep));
    iFFTs = [];
    for y = 1:length(yStep)
        for x = 1:length(xStep)

            subI = I(yStep(y):yStep(y)+sS-1,xStep(x):xStep(x)+sS-1);

            fftI = abs(fft(subI,[],1));
            mfI = mean(fftI,2);
            fftCut = mfI(2:fix(length(mfI)/2));
            iFFTs(size(iFFTs,1)+1,:) = fftCut;

    

        end
    end

    subplot(1,1,1)
    
    median(iFFTs,1)
    for i = 1:size(iFFTs,1)
       plot(iFFTs(i,:))
       hold on
    end
    
    medFFT = median(iFFTs,1);
    plot(medFFT,'r','LineWidth',2)
    ylim([0 max(medFFT) * 1.2])
    hold off
    
%% Get Criteria

startHigh = 1;
stopHigh = 20;
startLow = 25;
stopLow = size(medFFT,2);
    
%for i = 1 : size
bgFFT = mean(fftCut(startLow:stopLow));
fftMos(y,x) = sum(fftCut(startHigh:stopHigh)-bgFFT);

%     
%     subplot(2,1,1),
%     image(I)
%     subplot(2,1,2)
%     image(fftMos/100),
%     pause(.01)

    sortFM = sort(fftMos(:),'descend')/1000;
    grab = min(fix(length(sortFM) * minSigArea)+1,length(sortFM))-1;
    highFFT = median(sortFM(1:grab));
    medianFFT = median(sortFM);
    lowFFT = median(sortFM(length(sortFM)-grab:end));
    percentSat = (sum(I(:)==0)+ sum(I(:)==255))/(ys * xs) * 100;
    
%     
%     showFFT{i,1} = iNam{i};
%     showFFT{i,2} = highFFT(i);
%     showFFT{i,3} = medianFFT(i);
%     showFFT{i,4} = lowFFT(i);
%     showFFT{i,5} = percentSat(i);
% 
%     plot(sortFM),pause(.1)
%     hold on
    

% hold off
bar(highFFT)

%% find bad

bestHalf = sort(highFFT,'descend');
bestHalf = bestHalf(1:min(fix(length(bestHalf)/3)+1,length(bestHalf)));
targ  = mean(bestHalf);
bHstd = std(bestHalf);
cutOff = targ - (bHstd * 3);
badies = find(highFFT<cutOff);
badNames = {iNam{badies}}';
badFs = highFFT;
badFs(badFs>= cutOff) = 0;
bar(highFFT)
hold on
bar(badFs,'FaceColor','r','EdgeColor','r')
hold off

xlswrite([TPN 'badFocus.xls'],showFFT,'Results');
xlswrite([TPN 'badFocus.xls'],badNames,'BadPics');
save([TPN 'showFFT.mat'],'showFFT')

scatter(highFFT,percentSat)

profile off
