%%

%% Required input

clear all
colormap gray(256)

reps = 100;  % number of randomizations
reportOn = 1; % display confirmation of proper matrix after each randomization
imageOn = 0;
countMults = 1;
countCohorts = 0;

%% Network variables
prefRat = 0.8;  %proportion of picks in which synapses are chosen according to bias


%% get touch and synapse matrix
load('singlespinematrices2.mat')
%load('singlespinematrices_josh.mat')
synMat = singlespinesynapsematrix>0;
touchMat = singlespinetouchmatrix>0;

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

if countCohorts
    commonSyn = matchDend(synMat,dendMat);
    randCommonSyn = zeros(dendNum,dendNum,reps);
end


%% Make storage Matricies
sumSyn = double(synMat*0);

%% Run model
for r = 1:reps
    
    %% make networks
    if 0 
        makeSyn = randSteal(synMat,touchMat);
    elseif 0
        makeSyn = randStealWeightedAxons(synMat,touchMat);
    elseif 0
        makeSyn = randBiasSteal(synMat,touchMat,dendMat,prefRat);
        prefRats(r) = prefRat;
    elseif 0
        prefRat = 0;
        makeSyn = minOneSteal(synMat,touchMat,dendMat,prefRat);
    elseif 0
        makeSyn = pickFromBestTouch(synMat,touchMat);
    elseif 0
        makeSyn = randSortTouch(synMat,touchMat);
    elseif 0
        makeSyn = driftGraph(synMat,touchMat);
    elseif 1 
        makeSyn = driftSolverMultiTouch(synMat*0,synMat,touchMat,0);
    end
    
    %% need to run by picking from touche perspective
    
    
    %% Analyze
    if  countMults
        [randUniques(r) randMults(r)] = findUniques(makeSyn,dendMat);
    end
    
    if countCohorts
        randCommonSyn(:,:,r) = matchDend(makeSyn,dendMat);
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
    
    disp(sprintf('%02.1f%% finished',r/reps*100))
    
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

if countMults
    
    Puniques = sum((randUniques)<=realUniques)/reps;
    disp(sprintf('Total uniques P = %0.05f',Puniques))
    
    Pmults = sum((randMults)>=realMults)/reps;
    disp(sprintf('Total mults P = %0.05f',Pmults))
    hist(randMults,.5:1:max(randMults))
    hold on
    scatter(realMults,1,'r')
    hold off
    
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
    
    subplot(1,1,1)
    scatter(randMults,randAllCom,'b')
    hold on
    scatter(realMults,realAllCom,'r')
    hold off
    
    %% 
    showCross = zeros(max(round(randAllCom)),max(round(randMults)));
    crossInd = sub2ind(size(showCross),round(randAllCom),round(randMults));
    histInd = hist(crossInd, 0.5 : 1 : length(showCross(:)));
    showCross(:) = histInd;
%     showCross = repmat(showCross,[1 1 3]);
%     showCross(:,:,2) = showCross(:,:,2) * 5;
%         showCross(:,:,3) = showCross(:,:,3) *.2;
%         showCross(:,:,1) = showCross(:,:,1) * 0;
%     showCross(round(realAllCom),round(realMults),1) = 1000;
%     
    chowCross = showCross*150/max(showCross(:));
    cmap = jet(256);
    cmap(1,:) = 0;
    colormap(cmap)
    
    subplot(3,1,1)
    image(uint8(showCross))
    ylim([300,515])
    hold on
        title('redundant synapse number vs shared axon number')

    scatter(realMults,realAllCom,'r','linewidth',3)
    hold off
    
    subplot(3,1,2)
    bar(showCross(realAllCom,:))
    hold on
        title('redundant number of synapses given shared axon number')

    scatter(realMults,1,'r','linewidth',2 )
    hold off
    
    subplot(3,1,3)
    bar(showCross(:,realMults))
        title('Shared axons given real redundant number')

    hold on
    scatter(realAllCom,1,'r','linewidth',2)
    hold off
    
    %ylim([0 600])
    %%
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





