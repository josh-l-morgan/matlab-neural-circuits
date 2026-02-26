clear p* x*  beta

xTime = linspace(-500,0,50);
figure(1)
clf
subplot(1,2,1)
hold on
plot(xTime,spikes(i).lowContrastEarlySTAnormalized,'m')
plot(xTime,spikes(i).lowContrastLateSTAnormalized,'b')
plot(xTime,spikes(i).highContrastEarlySTAnormalized,'g')
plot(xTime,spikes(i).highContrastLateSTAnormalized,'r')
% legend('low contrast 0-5s','low contrast 20-30s',...
%'high contrast 0-5s','high contrast 20-30s') 
% legend('boxoff')
xlabel('delay (ms)')
ylabel('relative STA intensity')
axis tight

subplot(1,2,2)
hold on
plot(spikes(i).lowContrastLateInput2,spikes(i).lowContrastLateOutput*100,'b.')
plot(spikes(i).lowContrastEarlyInput2,spikes(i).lowContrastEarlyOutput*100,'m.')
plot(spikes(i).highContrastEarlyInput2,spikes(i).highContrastEarlyOutput*100,'g.')
plot(spikes(i).highContrastLateInput2,spikes(i).highContrastLateOutput*100,'r.')


subplot(1,2,2)
hold on
%highContrastEarly
M = max(spikes(i).highContrastEarlyOutput);
x = linspace(spikes(i).highContrastEarlyInput2(1),...
    spikes(i).highContrastEarlyInput2(end),1000);
beta = spikes(i).highContrastEarlyParameters;
p = normcdf ( (beta(1) * x + beta(2)), 0, 1);
y = M * p;
plot(x,y*100,'g')
% plot(x,y,'g')

%highContrastLate
x = linspace(spikes(i).highContrastLateInput2(1),...
    spikes(i).highContrastLateInput2(end),1000);
beta = spikes(i).highContrastLateParameters;
p = normcdf ( (beta(1) * x + beta(2)), 0, 1);
y = M * p;
plot(x,y*100,'r')
% plot(x,y,'r')

%lowContrastEarly
x = linspace(spikes(i).lowContrastEarlyInput2(1),...
    spikes(i).lowContrastEarlyInput2(end),1000);
beta = spikes(i).lowContrastEarlyParameters;
p = normcdf ( (beta(1) * x + beta(2)), 0, 1);
y = M * p;
plot(x,y*100,'m')

%lowContrastLate
x = linspace(spikes(i).lowContrastLateInput2(1),...
    spikes(i).lowContrastLateInput2(end),1000);
beta = spikes(i).lowContrastLateParameters;
p = normcdf ( (beta(1) * x + beta(2)), 0, 1);
y = M * p;
plot(x,y*100,'b')
% plot(x,y,'b.')
xlabel('Input')
ylabel('Output (Hz)')
axis tight

highDrive = (spikes(i).highContrastLateParameters(2)...
    /spikes(i).highContrastLateParameters(1))
lowEDrive = (spikes(i).lowContrastEarlyParameters(2)...
    /spikes(i).lowContrastEarlyParameters(1))
lowLDrive = (spikes(i).lowContrastLateParameters(2)...
    /spikes(i).lowContrastLateParameters(1))


highSens = spikes(i).highContrastLateParameters(1)
lowESens = spikes(i).lowContrastEarlyParameters(1)
lowLSens = spikes(i).lowContrastLateParameters(1)


