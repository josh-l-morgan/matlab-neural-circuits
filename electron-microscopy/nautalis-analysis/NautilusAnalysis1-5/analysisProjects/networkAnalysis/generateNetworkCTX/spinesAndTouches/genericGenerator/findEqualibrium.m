%%

%% Required input

clear all
colormap gray(256)

keepSynMat = 0; % 1 if start with real synaptic matrix
metaReps = 10000;
reps = 50;  % number of randomizations
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 0;
countMults = 0;
countMultTouch = 0;
countCohorts = 0;


%% Network variables
prefRat = 0.8;  %proportion of picks in which synapses are chosen according to bias


%% get touch and synapse matrix
 load('singlespinematrices2.mat')
%load('C:\Users\joshm\Documents\MATLAB\jlm_Code\models\stastiticalFarce\sF_results.mat')

if exist('spineparentids','var')
    synMat = singlespinesynapsematrix>0;
    touchMat = singlespinetouchmatrix>0;
    countMults = 1;
else
    synMat = sF_results.synMat;
    touchMat = sF_results.touchMat;
    spineparentids = 1:size(touchMat,2);
    countMultTouch = 1;
end
    

sortDend = sort(unique(spineparentids));
lookupDend(sortDend) = 1:length(sortDend);
parentDend = lookupDend(spineparentids);
dendMat = repmat(parentDend,[size(touchMat,1) 1]);


%% get per cell synapse numbers (constraints)
aSynNum = sum(synMat,2)';
sSynNum = sum(synMat,1);
histDend = [.5:max(dendMat(:))];
dendNum = max(dendMat(:));
synFrac = sum(synMat(:))/length(synMat(:));



%% Analyze real data

if countMults
    [realUniques realMults] = findUniques(synMat,dendMat);
end

if countMultTouch
    realMults = sum(synMat(:))-sum(synMat(:)>0);
    realUniques = sum(synMat(:)>0);
end

if countCohorts
    commonSyn = matchDend(synMat,dendMat);
    randCommonSyn = zeros(dendNum,dendNum,reps);
end


%% Make storage Matricies
sumSyn = double(synMat);


parfor mr = 1:metaReps
%% Run model
makeSyn = synMat * keepSynMat;
for r = 1:reps
    
    %% make networks
    makeSyn = driftSolverMultiTouch(makeSyn,synMat,touchMat,1);
    
    
    %% Analyze
    if  countMults
        [randUniques(mr,r) randMults(mr,r)] = findUniques(makeSyn,dendMat);
    end
    
    if countMultTouch
        randMults(mr,r) = sum(makeSyn(:))-sum(makeSyn(:)>0);
        randUniques(mr,r) = sum(makeSyn(:)>0);
    end
    
    if countCohorts
       [ sharedAx randCommonSyn(:,:,r)] = matchDend(makeSyn,dendMat);
    end
    
    
    %% Report
    if reportOn
        synDif = checkNewNetwork(synMat,makeSyn,touchMat);
        sumSyn = sumSyn + makeSyn;
    end
    
    %% Image
    if imageOn
        col = synMat*1000;
        col(:,:,2) = makeSyn * 1000;
        col(:,:,3) = touchMat*1000;
        image(uint8(col))
    end
    
    
end
    disp(sprintf('%02.1f%% finished',mr/metaReps*100))

end

%% Analyze distribution of synapses produced by model

onOld = sum(sum(sumSyn(synMat>0)));
onNew = sum(sum(sumSyn(synMat <1)));

newRat = onNew/(onOld+onNew)*100;
binSynReps = [0:.01:1];
synReps = (hist(sumSyn(touchMat>0)/reps,binSynReps))/sum(synMat(:));
bar(synReps)


%% Image Final
col = synMat*1000;
col(:,:,2) = sumSyn*10000;
col(:,:,3) = touchMat*1000;
image(uint8(col)),pause(.01)

%% Analyze Results

if countMults | countMultTouch
    
    Puniques = sum((randUniques)<=realUniques)/reps;
    disp(sprintf('Total uniques P = %0.05f',Puniques))
    
    Pmults = sum((randMults)>=realMults)/reps;
    disp(sprintf('Total mults P = %0.05f',Pmults))
    hist(randMults,.5:1:max(randMults))
    %plot(randMults')
    [sy sx] = find(randMults>-10,'.');
    
    hold on
    scatter(realMults,1,'r')
    hold off
    
    showMults = zeros(60,size(randMults,2));
    for i = 1:size(randMults,1)
        showInd = sub2ind(size(showMults),randMults(i,:)+1,1:size(randMults,2));
        showMults(showInd) = showMults(showInd) + 1;
    end
    showMults =showMults *3000/metaReps;
    image(showMults)
    hold on
    plot(mean(randMults,1),'b','LineWidth',5)
    hold on
    plot(ones(size(randMults,2),1)*(realMults+1),'r','LineWidth',5);
    hold off
    
    colMults = repmat(showMults,[1 1 3]);
    colMults(realMults,end,:) = colMults(realMults,end,:) + 100; 
        colMults(realMults,:,1) = colMults(realMults,:,1) + 100; 
subplot(1,1,1)
    %image(uint8(colMults))
    
    TPN = 'C:\Users\joshm\Documents\myWork\myProjects\NetworkAnalysisJournalClub\myNetworkGenerationPresentation\data\'
    imwrite(uint8(colMults),[TPN 'startFromReallJitLGN_1000.tif'])
    
    
    
    sortedUs = sort(randMults);
    threshMults = sortedUs(round(length(sortedUs)*.95));
    difMults = threshMults - realMults;
    percentChangeMults = difMults/realMults*100;
    meanRandMults = mean(randMults)
    
end

if countCohorts
    
    selfs = sub2ind(size(commonSyn),[1:dendNum],[1:dendNum]);
    for r = 1:reps
        grabCom = randCommonSyn(:,:,r);
        grabCom(selfs) = 0;
        randAllCom(r) = sum(grabCom(:))/2;
    end
    realCom = commonSyn;
    realCom(selfs) = 0;
    realAllCom = sum(realCom(:))/2;
    hist(randAllCom,[350:3:520])
    hold on
    scatter(realAllCom,0,'r')
    hold off
    
    Pcom = sum(randAllCom<=realAllCom)/reps
    
    meanRandCom = mean(randCommonSyn,3);
    meanRandCom(selfs) = 0;
    colCom = meanRandCom*40;
    colCom(:,:,2) = realCom * 40;
    colCom(:,:,3) = meanRandCom*80;
    image(uint8(colCom))
    hist(randAllCom,[351:5:520])
    
    
    % meanRandComSyn = mean(randCommonSyn,3);
    % for r = 1:reps
    %
    %     difRandCom(:,:,r) = randCommonSyn(:,:,r)-meanRandComSyn;
    %     %image(difRandCom(:,:,r)*100+100),pause(.1)
    % end
    % meanDifRandCom = mean(difRandCom,3);
    % randMoreCom = difRandCom * 0;
    % randMoreCom(difRandCom>0) = difRandCom(difRandCom>0);
    % randLessCom = difRandCom * 0;
    % randLessCom(difRandCom<0) = difRandCom(difRandCom<0);
    % difRealCom = commonSyn - meanRandComSyn;
    %
    %
    % realMoreCom = difRealCom(difRealCom>0);
    % realLessCom = difRealCom(difRealCom<0);
    %
    % randSumAllCom = squeeze(sum(sum(randCommonSyn,1),2));
    % randSumMoreCom = squeeze(sum(sum(randMoreCom,1),2));
    % randSumLessCom = squeeze(sum(sum(randLessCom,1),2));
    % realSumMoreCom = sum(sum(realMoreCom,1),2);
    % realSumLessCom = sum(sum(realLessCom,1),2);
    % realSumCommonSyn = sum(sum(commonSyn));
    %
    %
    % PlessCom = sum(randSumLessCom<=realSumLessCom)/reps
    % PMoreCom = sum(randSumMoreCom<=realSumMoreCom)/reps
    % PAllCom = sum(randSumAllCom<=realSumCommonSyn)/reps
    % hist(randSumAllCom,[1450.5:1:1700])
    
    
%     image(commonSyn*100)
%     clear comCol
%     comCol(:,:,1) = commonSyn;
%     comCol(:,:,2) = meanRandComSyn;
%     comCol(:,:,3) = commonSyn;
%     image(uint8(comCol*100))
    
end





