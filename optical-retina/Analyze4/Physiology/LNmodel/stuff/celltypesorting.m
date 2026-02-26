%% LOADING FILES AND GETTING NUMBER OF CHANNELS
clc
clear spikes
[filename,directoryname] = uigetfile;
load([directoryname filename])
cd('C:\Program Files\MATLAB\R2006a\work')
s = whos('spikes');
spikes.channel
nChannels = s.size(2)


%% CELL TYPE ASSIGNMENT
onBsFast = [];
onBtSlow =[];
offBsBiphasic = [];
offBtMonophasic =[];
led=[];
onOffDs =[];

if exist('ONBSFAST', 'var') && any(onBsFast)
    s = whos('ONBSFAST');
    nExperiments = s.size(2);
    ONBSFAST(nExperiments + 1).filename = filename;
    ONBSFAST(nExperiments + 1).spikes = spikes(onBsFast);
elseif any(onBsFast)
    ONBSFAST(1).filename = filename;
    ONBSFAST(1).spikes = spikes(onBsFast);
else
end

if exist('ONBTSLOW', 'var') && any(onBtSlow)
    s = whos('ONBTSLOW');
    nExperiments = s.size(2);
    ONBTSLOW(nExperiments + 1).filename = filename;
    ONBTSLOW(nExperiments + 1).spikes = spikes(onBtSlow);
elseif any(onBtSlow)
    ONBTSLOW(1).filename = filename;
    ONBTSLOW(1).spikes = spikes(onBtSlow);
else
end

if exist('OFFBSBIPHASIC', 'var') && any(offBsBiphasic)
    s = whos('OFFBSBIPHASIC');
    nExperiments = s.size(2);
    OFFBSBIPHASIC(nExperiments + 1).filename = filename;
    OFFBSBIPHASIC(nExperiments + 1).spikes = spikes(offBsBiphasic);
elseif any(offBsBiphasic)
    OFFBSBIPHASIC(1).filename = filename;
    OFFBSBIPHASIC(1).spikes = spikes(offBsBiphasic);
else
end

if exist('OFFBTMONOPHASIC', 'var') && any(offBtMonophasic)
    s = whos('OFFBTMONOPHASIC');
    nExperiments = s.size(2);
    OFFBTMONOPHASIC(nExperiments + 1).filename = filename;
    OFFBTMONOPHASIC(nExperiments + 1).spikes = spikes(offBtMonophasic);
elseif any(offBtMonophasic)
    OFFBTMONOPHASIC(1).filename = filename;
    OFFBTMONOPHASIC(1).spikes = spikes(offBtMonophasic);
else
end

if exist('LED', 'var') && any(led)
    s = whos('LED');
    nExperiments = s.size(2);
    LED(nExperiments + 1).filename = filename;
    LED(nExperiments + 1).spikes = spikes(led);
elseif any(led)
    LED(1).filename = filename;
    LED(1).spikes = spikes(led);
else
end

if exist('ONOFFDS', 'var') && any(onOffDs)
    s = whos('ONOFFDS');
    nExperiments = s.size(2);
    ONOFFDS(nExperiments + 1).filename = filename;
    ONOFFDS(nExperiments + 1).spikes = spikes(onOffDs);
elseif any(onOffDs)
    ONOFFDS(1).filename = filename;
    ONOFFDS(1).spikes = spikes(onOffDs);
else
end

