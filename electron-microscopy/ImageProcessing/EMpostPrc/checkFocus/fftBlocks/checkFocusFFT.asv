clear all
colormap gray(256)
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



%% Start Reading
for i = 1:min(10,length(iNam));
    i
I = double(imread([TPN iNam{i}]));
%fftI = abs(fftshift(fft2(I)));
fftI = abs(fft(I-mean(I(:)),[],1));
mfI = mean(fftI,2);
ffts(i,:) = mfI(2:fix(length(mfI)/2));


end


%% 
colU = ['b','r'];
for i = 1: size(ffts,1)
    useFFT = ffts(i,:);
    bgFFT = mean(useFFT(300:end));
    fftRat(i) = mean(useFFT(1:300))/bgFFT;
    fftTot(i) = sum(useFFT(1:20)-bgFFT);   
    if i == 1, c = 1; else c = 2; end
    plot(useFFT,colU(c)),i , pause
    hold on
 
end
hold off
plot(ffts(1,:)-ffts(2,:))
plot(mean(ffts,1))
plot((ffts(1,:)-ffts(2,:))./mean(ffts,1))



%bar(fftTot)
% 
% plot(fftTot/max(fftTot(:)))
% hold on
% plot(fftRat/max(fftRat(:)),'r')
% hold off


