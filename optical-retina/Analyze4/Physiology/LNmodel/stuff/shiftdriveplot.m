%% ONBSFAST
d=5;
unityLine = linspace(0,3,1000);
[m nExperiments] = size(ONBSFAST);
for k=1:nExperiments
    xSens = ONBSFAST(k).spikes.highContrastFit(:,1)./sensMed;
    ySens = ONBSFAST(k).spikes.lowContrastFit(:,1)./sensMed;
    xDrive = ONBSFAST(k).spikes.highContrastFit(:,2)./driveMed;
    yDrive = ONBSFAST(k).spikes.lowContrastFit(:,2)./driveMed;

    subplot(2,2,1)
    plot(unityLine, unityLine,'k')
    hold on
    plot(xSens, ySens, 'ro', 'MarkerSize', d)

    subplot(2,2,3)
    plot(unityLine, unityLine,'k')
    hold on
    plot(xDrive, yDrive, 'ro', 'MarkerSize', d)
end
clear x* y*

%% ONBTSLOW
[m nExperiments] = size(ONBTSLOW);
for k=1:nExperiments
    xSens = ONBTSLOW(k).spikes.highContrastFit(:,1)./sensMed;
    ySens = ONBTSLOW(k).spikes.lowContrastFit(:,1)./sensMed;
    xDrive = ONBTSLOW(k).spikes.highContrastFit(:,2)./driveMed;
    yDrive = ONBTSLOW(k).spikes.lowContrastFit(:,2)./driveMed;

    subplot(2,2,1)
    plot(xSens, ySens, 'rs', 'MarkerSize', d)

    subplot(2,2,3)
    plot(xDrive, yDrive, 'rs', 'MarkerSize', d)
end
clear x* y*

%% OFFBSBIPHASIC
[m nExperiments] = size(OFFBSBIPHASIC);
for k=1:nExperiments
    xSens = OFFBSBIPHASIC(k).spikes.highContrastFit(:,1)./sensMed;
    ySens = OFFBSBIPHASIC(k).spikes.lowContrastFit(:,1)./sensMed;
    xDrive = OFFBSBIPHASIC(k).spikes.highContrastFit(:,2)./driveMed;
    yDrive = OFFBSBIPHASIC(k).spikes.lowContrastFit(:,2)./driveMed;

    subplot(2,2,2)
    plot(unityLine, unityLine,'k')
    hold on
    plot(xSens, ySens, 'ro', 'MarkerSize', d)

    subplot(2,2,4)
    plot(unityLine, unityLine,'k')
    hold on
    plot(xDrive, yDrive, 'ro', 'MarkerSize', d)
end
clear x* y*

%% OFFBTMONOPHASIC
[m nExperiments] = size(OFFBTMONOPHASIC);
for k=1:nExperiments
    xSens = OFFBTMONOPHASIC(k).spikes.highContrastFit(:,1)./sensMed;
    ySens = OFFBTMONOPHASIC(k).spikes.lowContrastFit(:,1)./sensMed;
    xDrive = OFFBTMONOPHASIC(k).spikes.highContrastFit(:,2)./driveMed;
    yDrive = OFFBTMONOPHASIC(k).spikes.lowContrastFit(:,2)./driveMed;

    subplot(2,2,2)
    plot(xSens, ySens, 'rs', 'MarkerSize', d)

    subplot(2,2,4)
    plot(xDrive, yDrive, 'rs', 'MarkerSize', d)
end
clear x* y*
