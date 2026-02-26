%% Required input
% Every axon and spine must be able to make at least one synapse

clear all
colormap gray(256)

reps = 100;  % number of randomizations
displayInterval = 100000;  % number of iterations per randomaization before displaying status
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 1;

%% get touch and synapse matrix
load('singlespinematrices2.mat')
synMat = singlespinesynapsematrix>0;
touchMat = singlespinetouchmatrix>0;
spineparentids


%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);



%% Run model
tic
[ys xs] = size(touchMat);
synRat = sum(synMat(:))/sum(touchMat(:));
numTouch = sum(touchMat(:)>0);
touchInd = find(touchMat>0);
for r = 1:reps
    randTouch = rand(numTouch,1);
    randSyn = randTouch <= synRat;
    
    
    
end




