
%clear all
%generic connectivity skew test
%-joshm
tic
%% options
reps = 10000;   %how many repetitions to do

%% get data
'getting data'
% [preIn postIn] = parseBinDat;
cons = load('C:\Users\joshm\Documents\MyWork\myProjects\bobbysCTX\moreData_12+05+22\cons.mat')
preIn = cons.preIn;
postIn = cons.postIn;

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

multAx = find(hitNum>1 ); %Select population of axons to test
sampSize = length(multAx)

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
countSyn = 0;
    for a = 1:length(multAx) %run analysis for all axons making more than one synapse
        ax = multAx(a); %get multi axon id
        hitPost = post(axPos{ax}) %retrieve post partners for axon 
        countSyn = countSyn + length(hitPost);
        aHist = hist(hitPost,[.5:1:(postNum-.5)]); % convert to histogram of hits
        if max(aHist)>1 % if pre makes multiple synapses onto at least one post
            dDif =  aHist - postDist * hitNum(ax); % find difference betwen predicted and observed
             difD= fix(max(abs(dDif))); %record maximum difference over predicted (of the three dendrits)
        else
            difD = 0; % do difference if no convergence
        end
        
        observedDifs(a) = difD; % record skews
    end
subplot(2,1,1)
hist(observedDifs) %plot observed distribution of skews
[obHist xout] = hist(observedDifs);
obSkew = sum(observedDifs)/countSyn; %percent of displaced synapses

%% run
'running randomization'

allDifs = [];

aveSkew = zeros(reps,1);  %set up matrix to record average skews
difs = zeros(length(multAx),1);  %set up matrix to record skew for each axon with multiple synapses
aveHist = obHist*0;
for i = 1:reps
    %%Scramble posts
    %pos = rand(synNum,1); %generate list of random numbers
    %[a newPos] = sort(pos); %sort list to generate list of random positions
    %newDend = post(newPos); %use random positions to rearange posts
    newDend = post(fix(rand(synNum,1)*length(post))+1);
    for a = 1:length(multAx) %run each axon with more than one synapse
        ax = multAx(a); % retrive mult axon id
        hitPost = newDend(axPos{ax}); %find dendrits innervated by each axon
        aHist = hist(hitPost,[.5:1:(postNum-.5)]); %make histogram of partners
        if max(aHist)>1 % if multiple innervation or one or more dendrite
            dDif =  aHist - postDist * hitNum(ax); %measure skewness as before
            difD= fix(max(abs(dDif))); %take maximum
        else
            difD = 0; %no skew if no convergence
        end
        
        difs(a) = difD; %record diffrences for all multiple axons with more than one contact
    end
    allDifs = [allDifs difs];
    aveSkew(i) = sum(difs)/countSyn; %record averages for all runs
    aveHist = aveHist + hist(difs,xout);
    
    if mod(i,1000) == 0
        sprintf('ran %d of %d reps',i,reps) %report progress
    end
end

%% Analyze Results


P = mean(aveSkew >= obSkew);  %find P value
[num2str(length(multAx)) ' axons tested. Average skew from model was ' num2str(mean(aveSkew)) ...
    '. Observed skew was ' num2str(obSkew) '. P = ' num2str(P) ' for ' num2str(reps) ' reps.']

%% Survival curve

sortSkew = sort(aveSkew);
probs = .1:.001:1;
for i = 1:length(probs)
    disAtProb(i) = sortSkew(round( length( sortSkew )* probs(i) ))
end
percentObDisAtProb = (obSkew - disAtProb) * 100
percentObDisAt95 = (obSkew - sortSkew(round(length(sortSkew) * 0.975)))*100;

%% Display data

skewBins = [min(aveSkew):(max(aveSkew)-min(aveSkew))/100:max(aveSkew)]; 
subplot(3,1,1)
bar1=bar(xout, aveHist/(reps * preNum), 'FaceColor', 'b', 'EdgeColor', 'k'); 
set(bar1,'BarWidth',.8); 
hold on;
bar2=bar(xout, obHist/preNum, 'FaceColor', 'r', 'EdgeColor', 'k');
set(bar2,'BarWidth',.4); 
hold off; 
legend('randomized','observed')
XLabel('Skewness')
YLabel('proportion of sample')

subplot(3,1,2)
allHist = hist(aveSkew,skewBins)/reps;
bar4 = bar(skewBins,allHist,'b')
set(bar4,'BarWidth',length(allHist)/30); 
hold on
higherHist = hist(aveSkew(aveSkew>=obSkew),skewBins,'r')/reps;
bar3 = bar(skewBins,higherHist,'r')
set(bar3,'BarWidth',length(allHist)/30); 

scatter(obSkew,0,'r','LineWidth',3)
hold off
XLabel('Fraction of total synapses that are displaced')
YLabel('Proportion of randomizations')
legend('randomized','observed')


subplot(3,1,3)
plot(probs*100,percentObDisAtProb)
line([0 100],[0 0],'Color','k')
line([0 95],[percentObDisAt95 percentObDisAt95],'Color','r')
line([95 95],[0 percentObDisAt95],'Color','r')
XLabel('Percent Confidence')
YLabel('Percent of synapses displaced')
toc


%% Calc effect size.

meanObDif = mean(observedDifs);
effectSize = sum(allDifs(:) < meanObDif)/length(allDifs(:))
reverseES = sum(observedDifs>mean(allDifs(:)))/length(observedDifs)

%% Reverse P
pcrit = 0.05;




