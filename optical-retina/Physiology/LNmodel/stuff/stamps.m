% [fileName directoryName]=uigetfile;
% load(fileName)
% i=1;
clf

kernelTime = (-500:10);
autoTime = linspace(0,100,100);
diffuseTime = linspace(-0.5,8.5,180);
STA = spikes(i).STA;
STA(1,:)=0.5;
tpeak = find (std(STA) == max(std(STA)));
fpeak = mean(STA(:,tpeak-2:tpeak+2),2);
mg = mean(STA(:));
low = min(STA(:));
high = max(STA(:));



subplot(2,2,1)
plot(kernelTime, spikes(i).time.normkernel,'b')
xlabel('delay(ms)')
ylabel('relative STA intensity')
axis tight
title(fileName)



subplot(2,2,3)
plot(autoTime, spikes(i).autocorrelation, 'r')
xlabel('time(ms)')
ylabel('relative rate')
axis tight


subplot(2,2,4)
plot(diffuseTime, spikes(i).diffuse, 'k')
xlabel ('time(s)')
ylabel('firing rate (Hz)')
axis tight

subplot(2,2,2)
if (mg-low) >= (high-mg)
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([low, 1 - low]);
else
    snap = reshape(fpeak, 80, 60)';
    pcolor (snap)
    shading flat
    colormap gray
    caxis([1 - high, high]);
end
sigmastr = strcat('\sigma = ', num2str(spikes(i).space.sigma));
TTPstr = strcat('TTP = ', num2str(spikes(i).time.peak));
str = {sigmastr; TTPstr};
text(10,10, str)
axis off
title(spikes(i).channel)


