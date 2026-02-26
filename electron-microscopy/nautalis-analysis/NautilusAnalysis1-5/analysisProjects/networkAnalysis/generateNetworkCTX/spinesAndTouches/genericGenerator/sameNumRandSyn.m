%function[makeSyn] = sameNumRandSyn(synMat,TouchMat)

reps = 100000;

countMults = 1;
countCohorts = 1;


 %% Analyze
    tic
    if  countMults
        [reaUniques realMults] = findUniques(synMat,dendMat);
    end
    toc
    tic
    if countCohorts
        realAllCom = matchDend(synMat,dendMat);
    end
    toc
    




%% Syn Num
synNum = sum(synMat(:));
touchNum = sum(touchMat(:));

%% touch list

catTouch = [];
for i = 1:max(touchMat(:))
    catTouch = cat(1,catTouch,find(touchMat == i));
end
toc

%% make syns
%randTouch = rand(length(lastCat),1);
%synTouch = randTouch<=synProb;

%% make syns
randUniques = zeros(reps,1);
randMults = zeros(reps,1);
sharedAx = zeros(reps,1);
for r = 1 : reps
    disp(sprintf('finished %02.1f of %d',r,reps))
    tic
    makeSyn = synMat * 0;
    
    randTouch = randperm(touchNum);
    randSyn = catTouch(randTouch<= synNum);
    uniqueSyn = unique(randSyn);
    synVals = hist(randSyn,uniqueSyn);
    
    makeSyn(uniqueSyn) = synVals;
    toc
    
    
    %% Analyze
    tic
    if  countMults
        [randUniques(r) randMults(r)] = findUniques(makeSyn,dendMat);
    end
    toc
    tic
    if countCohorts
        sharedAx(r) = matchDend(makeSyn,dendMat);
    end
    toc
    
end

%%  %% 
    showCross = zeros(max(round(sharedAx)),max(round(randMults)));
    crossInd = sub2ind(size(showCross),round(sharedAx),round(randMults));
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
    ylim([200,515])
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
    
    %