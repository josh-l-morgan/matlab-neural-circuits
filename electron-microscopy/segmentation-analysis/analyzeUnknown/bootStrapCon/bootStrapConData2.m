
%clear all
%generic connectivity skew test
%-joshm
tic
%% options
reps = 10000;   %how many repetitions to do

%% get data
'getting data'
[preIn postIn] = parseBinDat;

%% make dummy data to test model
% preIn = fix(rand(100,1)*90) + 1;
% postIn = fix(rand(100,1)*3)+1;

%% renumber IDs to insure every integer is used up to the number of pre and post cells

ids = sort(unique(preIn));  %get list of pre Ids
newIds = 1:length(ids); % make new lest of all integer IDs
lookupIds(ids) = newIds; % position new Ids into lookup table using old IDs as positions
pre = lookupIds(preIn); % use lookup tabel to transorm previous Ids into new Ids

ids = sort(unique(postIn)); %repeat of post partners
newIds = 1:length(ids);
lookupIds(ids) = newIds;
post = lookupIds(postIn);

%% get data information
preNum = length(unique(pre));  %number of presynaptic cells
postNum = length(unique(post)); %number of post synaptic cells
synNum = length(pre);  %number of synapses

postSites = hist(post,[.5:1:max(post)-.5]); %find number of sites on each post synaptic cell
postDist = postSites/synNum; %find ratio of sites given to each post cell
hitNum = hist(pre,[.5:1:max(pre)-.5]);  %find number of sites on each presynaptic cell
preDist = hitNum/synNum; % find ration of sites given to each pre cell. 

multAx = find(hitNum>1); %Identify axons making multiple synapses


%% map positions of all axons into cell array

for a = 1:preNum;
    axPos{a} = find(pre == a); %find each axon and stick its positions into cell array
end

%% measure observed skews
%%Skew is measured by first finding the expected innervation of post
%%synaptic cells based on the number of synapses for that axon and the
%%number of post synaptic sites.  This prediction is then subtracted from
%%the observed distribution.  The maximum difference (ie for the preferentially
%%innervated dendrite) is recorded as the axon's skew (difD).  
'measureing observed skews'
observedDifs = zeros(length(multAx),1); %create matix of observations
    for a = 1:length(multAx) %run analysis for all axons making more than one synapse
        ax = multAx(a); %get multi axon id
        hitPost = post(axPos{ax}); %retrieve post partners for axon 
        aHist = hist(hitPost,[.5:1:(postNum-.5)]); % convert to histogram of hits
        if max(aHist)>1 % if pre makes multiple synapses onto at least one post
            dDif =  aHist - postDist * hitNum(ax); % find difference betwen predicted and observed
            difD= max(dDif); %record maximum difference over predicted (of the three dendrits)
        else
            difD = 0; % do difference if no convergence
        end
        
        observedDifs(a) = difD; % record skews
    end
subplot(2,1,1)
hist(observedDifs) %plot observed distribution of skews
obSkew = mean(observedDifs); %average skew

%% run
'running randomization'

aveSkew = zeros(reps,1);  %set up matrix to record average skews
difs = zeros(length(multAx),1);  %set up matrix to record skew for each axon with multiple synapses
for i = 1:reps
    %%Scramble posts
    pos = rand(synNum,1); %generate list of random numbers
    [a newPos] = sort(pos); %sort list to generate list of random positions
    newDend = post(newPos); %use random positions to rearange posts
    
    for a = 1:length(multAx) %run each axon with more than one synapse
        ax = multAx(a); % retrive mult axon id
        hitPost = newDend(axPos{ax}); %find dendrits innervated by each axon
        aHist = hist(hitPost,[.5:1:(postNum-.5)]); %make histogram of partners
        if max(aHist)>1 % if multiple innervation or one or more dendrite
            dDif =  aHist - postDist * hitNum(ax); %measure skewness as before
            difD= max(dDif); %take maximum
        else
            difD = 0; %no skew if no convergence
        end
        
        difs(a) = difD; %record diffrences for all multiple axons with more than one contact
    end
    aveSkew(i) = mean(difs); %record averages for all runs
    
    if mod(i,1000) == 0
        sprintf('ran %d of %d reps',i,reps) %report progress
    end
end

%% Analyze Results

%%Display histogram of montecarlo results
skewBins = [min(aveSkew):(max(aveSkew)-min(aveSkew))/100:max(aveSkew)]; 
subplot(2,1,2)
hist(aveSkew,skewBins);

P = mean(aveSkew > obSkew);  %find P value
['Average skew from model was ' num2str(mean(aveSkew)) ...
    '. Observed skew was ' num2str(obSkew) '. P = ' num2str(P) ' for ' num2str(reps) ' reps.']



toc
