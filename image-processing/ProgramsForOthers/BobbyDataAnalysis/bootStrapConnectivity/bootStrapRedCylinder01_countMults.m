
%clear all
%generic connectivity skew test
%-joshm
tic
%% options
 load('C:\Users\joshm\Documents\MATLAB\models\bootStrapConnectivity\redCylinderExcite.mat')


reps = 10000;   %how many repetitions to do

%% get data
'getting data'
% [preIn postIn] = parseBinDat;

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

%%  count real mults
histX = [0:1:10];
maxAx = max(preIn);
offSetDend = maxAx+1;
shiftPostIn = postIn * offSetDend;
combSyn = shiftPostIn+preIn;
realCountAx = hist(combSyn,1:1:max(combSyn));
subCountAx = realCountAx(realCountAx>0);
realHistAx = hist(subCountAx,histX);
bar(histX,realHistAx); pause(.01)
%% run
'running randomization'

allDifs = [];

aveSkew = zeros(reps,1);  %set up matrix to record average skews
difs = zeros(length(multAx),1);  %set up matrix to record skew for each axon with multiple synapses

for i = 1:reps
    %%Scramble posts
    
%        pos = rand(synNum,1);
%     [a newPos] = sort(pos);
%     newDend = post(newPos);
    
    
    newDend = post(fix(rand(synNum,1)*length(post))+1);
    newAx = pre(fix(rand(synNum,1)*length(pre))+1);
    
    
    %newDend = post(fix(rand(synNum,1)*length(post))+1);
    %%  count real mults
    shiftPostIn = newDend * offSetDend;
    combSyn = shiftPostIn+preIn;
    CountAx = hist(combSyn,1:1:max(combSyn));
    subCountAx = CountAx(CountAx>0);
    HistAx(i,:) = hist(subCountAx,histX);
   

    %bar(histX,HistAx(i,:)); pause(.01)
    if mod(i,1000) == 0
        sprintf('ran %d of %d reps',i,reps) %report progress
    end
end

%% Analyze Results
actual(1,1) = realHistAx(3);
actual(2,1) = realHistAx(4);
actual(3,1) = sum(realHistAx(5:end));

testRes(1,:)= HistAx(:,3);
testRes(2,:) = HistAx(:,4);
testRes(3,:)= sum(HistAx(:,5:end),2);

sortRes = sort(testRes,2,'ascend');

Pthresh = 0.05;
testNum = 3;
acceptAt = Pthresh/testNum;

clear P t95
titles={'TWOS','THREES','MORE THAN THREE'};
for i = 1:3
   P(i,1) = sum(sortRes(i,:) >= actual(i))/reps; 
   cumulativePass(i,:) = sortRes(i,:) < actual(i);
   t95(i,1) = sortRes(i,round(reps*(1-acceptAt/testNum)));
   subplot(3,2,i*2-1)
   resHistX = [0:1:(actual(i)+1) * 2];

   histRes = hist(sortRes(i,:),resHistX)/reps;
   bar(resHistX,histRes)
   hold on
   scatter(t95(i),0,'r')
   scatter(actual(i),0,'g','*')
   title(titles{i})
   legend('model',sprintf('%2.1f%% threshold',(1-acceptAt)*100),...
       sprintf('actual (P = %1.5f)',P(i)))
   
        XLabel('number of axons with given convergence')
    YLabel('proportion of ramdomizations')
   hold off
end



%% Plot distributions
meanHistAx = mean(HistAx,1);
subplot(3,2,[2 4 6])
bar(histX,meanHistAx)
hold on
bar(histX,realHistAx,.3,'FaceColor','g')
% 
% U = meanHistAx*0;
% U(3:4) = t95(1:2)-meanHistAx(1:2)';
% errorbar(histX,meanHistAx,meanHistAx,U);

legend('model distribution', 'actual distribution')
XLabel('number of axons with given convergence')
YLabel('average / actual number of axons')
hold off

somePass = sum(cumulativePass,1)>0;
cumulativeP = sum(somePass)/reps

P
actual
t95


