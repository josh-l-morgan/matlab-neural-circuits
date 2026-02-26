clear p* x* high* low*
k=1;
spikes = OFFBTMONOPHASIC(k).spikes;
xTime = linspace(-500,0,50);

clf


subplot(1,2,1)
plot(xTime,spikes(i).lowContrastEarlySTAnormalized,'m')
hold on
plot(xTime,spikes(i).lowContrastLateSTAnormalized,'b')
plot(xTime,spikes(i).highContrastEarlySTAnormalized,'g')
plot(xTime,spikes(i).highContrastLateSTAnormalized,'r')
% legend('low contrast 0-5s','low contrast 20-30s','high contrast 0-5s','high contrast 20-30s') 
% legend('boxoff')
xlabel('delay (ms)')
ylabel('relative STA intensity')
axis tight



subplot(1,2,2)
plot(spikes(i).lowContrastLateInput,spikes(i).lowContrastLateOutput*100,'b*')
hold on
plot(spikes(i).lowContrastEarlyInput,spikes(i).lowContrastEarlyOutput*100,'m*')
plot(spikes(i).highContrastEarlyInput,spikes(i).highContrastEarlyOutput*100,'g*')
plot(spikes(i).highContrastLateInput,spikes(i).highContrastLateOutput*100,'r*')

%highContrastLate
xhighContrastLate = linspace(spikes(i).highContrastLateInput(1),spikes(i).highContrastLateInput(end),1000);
p1 = normcdf (xhighContrastLate,spikes(i).highContrastLateParameters(2), spikes(i).highContrastLateParameters(1));
highContrastLateFit =  p1 * spikes(i).highContrastLateParameters(3);
plot(xhighContrastLate,highContrastLateFit*100,'r')

%highContrastEarly
xhighContrastEarly = linspace(spikes(i).highContrastEarlyInput(1),spikes(i).highContrastEarlyInput(end),1000);
p2 = normcdf (xhighContrastEarly * spikes(i).highContrastEarlyParameters(1) + spikes(i).highContrastEarlyParameters(2),...
    spikes(i).highContrastLateParameters(2), spikes(i).highContrastLateParameters(1));
highContrastEarlyFit = p2 * spikes(i).highContrastLateParameters(3);
plot(xhighContrastEarly,highContrastEarlyFit*100,'g')

%lowContrastLate
xlowContrastLate = linspace(spikes(i).lowContrastLateInput(1),spikes(i).lowContrastLateInput(end),1000);
p3 = normcdf (xlowContrastLate * spikes(i).lowContrastLateParameters(1) + spikes(i).lowContrastLateParameters(2),...
    spikes(i).highContrastLateParameters(2), spikes(i).highContrastLateParameters(1));
lowContrastLateFit = p3 * spikes(i).highContrastLateParameters(3);
plot(xlowContrastLate,lowContrastLateFit*100,'b')

%lowContrastEarly
xlowContrastEarly = linspace(spikes(i).lowContrastEarlyInput(1),spikes(i).lowContrastEarlyInput(end),1000);
p4 = normcdf (xlowContrastEarly * spikes(i).lowContrastEarlyParameters(1) + spikes(i).lowContrastEarlyParameters(2),...
    spikes(i).highContrastLateParameters(2), spikes(i).highContrastLateParameters(1));
lowContrastEarlyFit = p4 * spikes(i).highContrastLateParameters(3);
plot(xlowContrastEarly,lowContrastEarlyFit*100,'m')

xlabel('Input')
ylabel('Output (Hz)')
axis tight