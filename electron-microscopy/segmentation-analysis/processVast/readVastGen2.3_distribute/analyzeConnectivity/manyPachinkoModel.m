
clear all
MPN = 'D:\LGNs1\Segmentation\VAST\S8\joshm\matObjects\matOut_14+05+28\'
TPN = [MPN 'manyTerSubs\'];
reps = 100000;


downSamps = [ 1 1 1 1 1  2 4 4  8 8 8  16 32 64]
lookDists = [ 0 1 2 3 4  4 2 4  4 8 16 16 16 16]

downSamps = [ 1  1 4 8  64]
lookDists = [ 0  4 4  8  16]


downSamps = [  4 ]
lookDists = [  4 ]

lookMicrons =  .03 * 4  + lookDists .* downSamps * 8 * 0.03

subplot(1,1,1)
for s = 1: length(lookDists)
    
downSamp = downSamps(s);
lookDist = lookDists(s);
fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)

load([TPN fileName]);
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;

%% Make overlap lookup table

sumOverlap = sum(touchMat(:));
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



subplot(length(lookDists),1,length(lookDists)-s+1)

if 0
meanRandSynDist = mean(randSynDist,1);
bar([meanRandSynDist; realSynDist']')
text(200,10,sprintf('%.2f micron dilation',lookMicrons(s)))

pause(.01)
else

%% 
divHist = (0:.02:3);
histRandSqrC = hist(randSqrC,divHist);
bar(divHist,histRandSqrC/max(histRandSqrC))
hold on
scatter(realSqrC,0.01,'r','lineWidth',10);
hold off
xlim([min(divHist) max(divHist)])
text(.3,.5,sprintf('%.2f micron dilation',lookMicrons(s)))

pause(.01)
end


manyRandSqrC(s).randSqrC = randSqrC;
manyRandSqrC(s).realSqrC = realSqrC;
manyRandSqrC(s).lookDist = lookDist;
manyRandSqrC(s).downSamp = downSamp;




end

%save([TPN 'manyPachinko.mat'],'manyRandSqrC')
%%

for i = 1:length(manyRandSqrC)
    lookMicrons(i) =  .03 * 4  + manyRandSqrC(i).lookDist * manyRandSqrC(i).downSamp * 8 * 0.03;
    meanRandSqrC(i) = mean(manyRandSqrC(i).randSqrC);
end
% 
% plot(lookMicrons,meanRandSqrC)
% %bar(lookMicrons,meanRandSqrC)
% xlim([-1 100])
% hold on
% scatter(0,manyRandSqrC(1).realSqrC,'r')
% hold off
