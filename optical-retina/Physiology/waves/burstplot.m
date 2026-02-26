function burstplot(waves, episode, channels)
clf
nChannels = length(waves);
% nChannels = 5;
hAxes = axes;


for i = 1:nChannels
    spikeTrain = waves(i).spikeTrain;
    spikeTrain(spikeTrain < episode(1)) = [];
    spikeTrain(spikeTrain > episode(2)) = [];
    nSpikes = length(spikeTrain);
    burstTime = waves(i).burstData(:,1:2);
    [nBursts nCols] = size(burstTime);
    waveTime = waves(i).waveData(:,1:2);
    [nWaves nCols] = size(waveTime);
    for j = 1:nSpikes
        plot(hAxes, [spikeTrain(j) spikeTrain(j)], [i i+0.2],'k');
        hold on
    end
    for j = 1:nBursts
        plot(hAxes, [burstTime(j,1) burstTime(j,2)], [i-0.2 i-0.2],'k')
        hold on
    end
    for j = 1:nWaves
        plot(hAxes, [waveTime(j,1) waveTime(j,2)], [i-0.5 i-0.5],'r')
        hold on
    end
end
set(    hAxes,...
    'YTickMode', 'manual',...
    'YTickLabelMode', 'manual',...
    'YTick', (1:nChannels),...
    'YTickLabel',channels');

axis([episode(1) episode(2) 0 nChannels+1]);

