%function[] = springModel(uselist,modPar);

% %%
% F = -kX
% k = stiffness
% X = proportional to distance;

%% Get data
seedList = [108 201];
useList = obI2cellList_seedInput(obI,seedList);
seedPref = seedPreferences(seedList,useList);

%% set variables
k = 1;
noise = 1;
damp = .9;
repulse = 10;
reps = 1000;
fsize = 100;

%% set fixed
springDir = 'D:\LGNs1\Analysis\springDat\'
wind = [0 fsize; 0 fsize];

%% configure nodes

postIDs = useList.postList;
preIDs = useList.preList;
postNum = length(postIDs);
preNum = length(preIDs);

postSynPref = seedPref.sharedSynNorm(1,:)./sum(seedPref.sharedSynNorm,1);
preSynPref = seedPref.ax2seed(1,:)./sum(seedPref.ax2seed,1);

con = useList.con;

postX = rand(postNum,1)*100;
postY = rand(postNum,1)*100;
postXV = zeros(postNum,1);
postYV = zeros(postNum,1);


preX = rand(preNum,1)*fsize;
preY = rand(preNum,1)*fsize;
preXV = zeros(preNum,1);
preYV = zeros(preNum,1);



%% run springs

for r = 1:reps
    disp(sprintf('running %d of %d',r,reps))
    
    
    preXmat = repmat(preX,[1,postNum]);
    preYmat = repmat(preY,[1,postNum]);
    postXmat = repmat(postX',[preNum,1]);
    postYmat = repmat(postY',[preNum,1]);

    %%Independent dim calc
    difX = preXmat - postXmat;
    difY = preYmat - postYmat;
    
    forceX = difX .* con * k;
    forceY = difY .* con * k;
    
    
    
    
    preForceX = sum(forceX,2);
    postForceX = sum(forceX,1);
    preForceY = sum(forceY,2);
    postForceY = sum(forceY,1);
    
    postXV = (postXV + postForceX') * damp;
    postYV = (postYV + postForceY'* damp);
    preXV = (preXV - preForceX)* damp;
    preYV = (preYV - preForceY)* damp;
    
    
    postX = postX + postXV;
    postY = postY + postYV;
    preX = preX + preXV;
    preY = preY + preYV;
    
    
    
    scatter(preX,preY,'k','.')
    hold on
    scatter(postX,postY,'k','o')
    hold off
    
    ylim([wind(1,1) wind(1,2)])
    xlim([wind(2,1) wind(2,2)])
    
    
    pause(.1)
    
end





















