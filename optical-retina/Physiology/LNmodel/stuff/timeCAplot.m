%% ONBSFAST
d=5;
peakLine = linspace(0,300,1000);
deltaLine = linspace (0,200,1000);
nExperiments = size(ONBSFAST,2);
for k=1:nExperiments
    xPeak = -ONBSFAST(k).spikes.highContrastTimePoint(:,2);
    yPeak = -ONBSFAST(k).spikes.lowContrastTimePoint(:,2);
    xDelta = ONBSFAST(k).spikes.highContrastTimePoint(:,3);
    yDelta = ONBSFAST(k).spikes.lowContrastTimePoint(:,3);

    subplot(2,2,1)
    plot(peakLine, peakLine,'k')
    hold on
    plot(xPeak, yPeak, 'k^', 'MarkerSize', d)

    subplot(2,2,3)
    plot(deltaLine, deltaLine,'k')
    hold on
    plot(xDelta, yDelta, 'k^', 'MarkerSize', d)
end
clear x* y*

%% ONBTSLOW
nExperiments = size(ONBTSLOW,2);
for k=1:nExperiments
    xPeak = -ONBTSLOW(k).spikes.highContrastTimePoint(:,2);
    yPeak = -ONBTSLOW(k).spikes.lowContrastTimePoint(:,2);
    xDelta = ONBTSLOW(k).spikes.highContrastTimePoint(:,3);
    yDelta = ONBTSLOW(k).spikes.lowContrastTimePoint(:,3);

    subplot(2,2,1)
    plot(xPeak, yPeak, 'ks', 'MarkerSize', d)

    subplot(2,2,3)
    plot(xDelta, yDelta, 'ks', 'MarkerSize', d)
end
clear x* y*

%% OFFBSBIPHASIC
nExperiments = size(OFFBSBIPHASIC,2);
for k=1:nExperiments
    xPeak = -OFFBSBIPHASIC(k).spikes.highContrastTimePoint(:,2);
    yPeak = -OFFBSBIPHASIC(k).spikes.lowContrastTimePoint(:,2);
    xDelta = OFFBSBIPHASIC(k).spikes.highContrastTimePoint(:,3);
    yDelta = OFFBSBIPHASIC(k).spikes.lowContrastTimePoint(:,3);

    subplot(2,2,2)
    plot(peakLine, peakLine,'k')
    hold on
    plot(xPeak, yPeak, 'k^', 'MarkerSize', d)

    subplot(2,2,4)
    plot(deltaLine, deltaLine,'k')
    hold on
    plot(xDelta, yDelta, 'k^', 'MarkerSize', d)
end
clear x* y*

%% OFFBTMONOPHASIC
nExperiments = size(OFFBTMONOPHASIC,2);
for k=1:nExperiments
    xPeak = -OFFBTMONOPHASIC(k).spikes.highContrastTimePoint(:,2);
    yPeak = -OFFBTMONOPHASIC(k).spikes.lowContrastTimePoint(:,2);
    xDelta = OFFBTMONOPHASIC(k).spikes.highContrastTimePoint(:,3);
    yDelta = OFFBTMONOPHASIC(k).spikes.lowContrastTimePoint(:,3);

    subplot(2,2,2)
    plot(xPeak, yPeak, 'ks', 'MarkerSize', d)

    subplot(2,2,4)
    plot(xDelta, yDelta, 'ks', 'MarkerSize', d)
end
clear x* y*
