
%clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];
load([TPN 'shiftTer_Ds8_Ds1_Look0_Dim2.mat']);

shiftVec = shiftTer.shiftVec;

reps = 10000;

lookDist = shiftTer.lookDist;
downSamp = shiftTer.downSamp;
lookMicrons =  .03 * 4  + lookDist .* downSamp * 8 * 0.03

subplot(1,1,1)
for s = 1: length(shiftVec)
   disp(sprintf('shift %d of %d',s,length(shiftVec))) 
%     
% downSamp = downSamps(s);
% lookDist = lookDists(s);
% fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)
% 
% load([TPN fileName]);
synMat = shiftTer.synMat;
touchMat = shiftTer.shiftMat(:,:,s);

%% Make overlap lookup table

sumOverlap = sum(touchMat(:));
ciOverlap(s) = squaredCluster(touchMat(:));


%{

[matY matX] = size(touchMat);
linMat = touchMat(:);
sumOverlap = sum(linMat);
lookupPair = zeros(sumOverlap,1);
prev = 0;
for i = 1:length(linMat);
    toEnd = prev + linMat(i);
    lookupPair(prev+1:toEnd) = i;
    prev = toEnd;
end
    
%% Make linear prediction
numSyn = sum(synMat(:));
predMat = touchMat * numSyn/sumOverlap;

realSqrC = squaredCluster(synMat);
realSynDist = sort(synMat(:),'descend');


%% Run model

randSynDist = zeros(reps,length(synMat(:)));
for r = 1: reps

    
makeSyn = synMat * 0;
pickRand = ceil(rand(numSyn,1)*sumOverlap);
pickedTouch = lookupPair(pickRand);
histSyn = hist(pickedTouch,[1:1:length(linMat)]);
makeSyn(:) = histSyn;

randSqrC(r) = squaredCluster(makeSyn);
randSynDist(r,:) = sort(makeSyn(:),'descend');

%image(makeSyn*10),pause(.01)

end

overlapClust(s) = 

subplot(length(shiftVec),1,length(shiftVec)-s+1)

if 0
meanRandSynDist = mean(randSynDist,1);
bar([meanRandSynDist; realSynDist']')
text(200,10,sprintf('%.2f micron dilation',lookMicrons(s)))

pause(.01)
else

%% 
divHist = (0:.05:11);
histRandSqrC = hist(randSqrC,divHist);
bar(divHist,histRandSqrC/max(histRandSqrC))
hold on
scatter(realSqrC,0.01,'r','lineWidth',10);
hold off
xlim([min(divHist) max(divHist)])
text(.3,.5,sprintf('%.2f micron dilation',shiftVec(s)))

pause(.01)
end

%}

manyRandSqrC(s).randSqrC = randSqrC;
manyRandSqrC(s).realSqrC = realSqrC;
manyRandSqrC(s).shiftVec = shiftVec;
manyRandSqrC(s).downSamp = downSamp;




end

ciOverlap
%save([TPN 'manyPachinko.mat'],'manyRandSqrC')
%%
% 
% for i = 1:length(manyRandSqrC)
%     %lookMicrons(i) =  .03 * 4  + manyRandSqrC(i).lookDist * manyRandSqrC(i).downSamp * 8 * 0.03;
%     
%     meanRandSqrC(i) = mean(manyRandSqrC(i).randSqrC);
% end

subplot(1,1,1)

shiftUm = shiftTer.shiftUm;

plot(shiftTer.shiftUm,ciOverlap,'LineWidth',3)
%bar(shiftVec,meanRandSqrC)
xlim([min(shiftUm) max(shiftUm)])
ylim([0 3])
hold on
scatter(0,manyRandSqrC(1).realSqrC,'r','LineWidth',3)
xlabel('Micrometer shift in axon position')
ylabel('Cluster index')
hold off
