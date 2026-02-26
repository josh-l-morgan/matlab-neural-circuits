%% Analyze FFT of subregions to check focus
%%Currently doing FFT only in Y dimension

clear all
colormap gray(256)

sS = 101;  % Define pixel length of area used for FFT sampling
minSigArea = .1; % proportion of image that must have signal

'Get good image'
[TFN TPN] =GetMyFile('Get Good Image');
riNam{1} = [TPN TFN];
'Get marginal image'
[TFN TPN] = GetMyFile('Get Marginal Image');
riNam{2} = [TPN TFN];


TPN = GetMyDir;


dTPN = dir(TPN);
dTPN = dTPN(3:end);

iNam = [];
for i = 1:length(dTPN)
    nam = dTPN(i).name;
    if length(nam)>4
        if strcmp(nam(end-3:end), '.tif' )
            iNam{length(iNam)+1} = nam;
        end
    end
end


%% Read Reference images
hold off
for i = 1:length(riNam);
    sprintf('running reference %d of %d',i,length(riNam))
    Ir = double(imread(riNam{i}));
    [mfIs, listYX] =funSubFFT(Ir);
    allMFIs(i,:) = mean(mfIs,1);
end
subplot(2,1,1)
for i = 1:size(allMFIs,1)
   plot(allMFIs(i,:))
   hold on
end
hold off

difMFI = ((allMFIs(1,:) - allMFIs(2,:)));

subplot(2,1,2)
plot(difMFI)
pause(.1)

aboveSTD = find(difMFI> (mean(difMFI(end-20:end))+std(difMFI(end-20:end))*3));
peakDif = find(difMFI == max(difMFI),1);
useMFI = aboveSTD(aboveSTD>peakDif);

goodMFI = allMFIs(1,useMFI);

%% Read data images
for i = 1:length(iNam);
    sprintf('running image %d of %d',i,length(iNam))
    Ir = double(imread([TPN iNam{i}]));
    [mfIs, listYX] =funSubFFT(Ir);
    %allMFIs(i,:) = mean(mfIs,1);
    
    
     fftCut = mfIs(:,useMFI);
     for m = 1:size(fftCut,1)
        fftM(m) = mean(fftCut(m,:)-goodMFI);
     end
     mys = max(listYX(:,1)); mxs = max(listYX(:,2));
     Iind = sub2ind([mys,mxs],listYX(:,1),listYX(:,2));
     fftMos = zeros(mys,mxs);
     fftMos(Iind) = fftM;
    
     subplot(2,1,1)
     image(Ir)
     subplot(2,1,2)
     image(fftMos*300),pause(.1)
        
     
     
    sortFM = sort(fftMos(:),'descend')/1000;
    grab = min(fix(length(sortFM) * minSigArea)+1,length(sortFM))-1;
    highFFT(i) = median(sortFM(1:grab));
    medianFFT(i) = median(sortFM);
    lowFFT(i) = median(sortFM(length(sortFM)-grab:end));
    percentSat(i) = (sum(Ir(:)==0)+ sum(Ir(:)==255))/(size(Ir,1) * size(Ir,2)) * 100;
        
    showFFT{i,1} = iNam{i};
    showFFT{i,2} = highFFT(i);
    showFFT{i,3} = medianFFT(i);
    showFFT{i,4} = lowFFT(i);
    showFFT{i,5} = percentSat(i);
% 
%     plot(sortFM),pause(.1)
%     hold on
     
end

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

%% write Exel
save([TPN 'showFFT.mat'],'showFFT')

if isempty(badNames),badNames = {'No failures found'};end
xlswrite([TPN 'badFocus.xls'],showFFT,'Results');
xlswrite([TPN 'badFocus.xls'],badNames,'BadPics');

%%Sort by highFFT
sortFFT = showFFT;
[sorted grab] = sort(highFFT,'ascend');
for i = 1: size(showFFT,1)
    sortFFT(i,:) = showFFT(grab(i),:);
end
xlswrite([TPN 'badFocus.xls'],sortFFT,'sortedResults');
%scatter(highFFT,percentSat)

