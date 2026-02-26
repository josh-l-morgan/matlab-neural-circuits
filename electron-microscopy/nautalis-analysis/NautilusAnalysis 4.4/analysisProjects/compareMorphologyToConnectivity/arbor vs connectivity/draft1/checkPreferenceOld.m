function[prefStat] = checkPreference(sourceCon,predSyn,useCon,seedTarg,randomizePick)


if randomizePick
    pickPref = randperm(size(sourceCon,1)); %choose new order of axons for determining seed preference
else
   pickPref = 1:size(sourceCon,1); 
end



%% find ax prefs without exclusions
for s = 1:length(seedTarg)
    prePrefCon(:,:,s) = repmat(sourceCon(pickPref,seedTarg(s)),[1 size(sourceCon,2)]);
%    subplot(2,1,s)
%     image(prePrefCon(:,:,s)*10)
end

%% Calculate what the preference of each post would be, absent a particular pre

postPrefMinusPre = zeros(size(sourceCon,1),size(sourceCon,2),2);
for i = 1: size(sourceCon,1)
    usePre = setdiff(1:size(sourceCon,1),i);
    usePick = pickPref(usePre);
    sampCon = sourceCon(usePre,:);
    prefCon = sourceCon(usePick,:);
    for s = 1:length(seedTarg)
        simMat = sqrt(sampCon .* repmat(prefCon(:,seedTarg(s)),[1 size(sourceCon,2)]));
        postPrefMinusPre(i,:,s) = sum(simMat,1);
    end
end
% 
% subplot(2,2,1)
% image(postPrefMinusPre(:,:,1));
% subplot(2,2,2)
% image(postPrefMinusPre(:,:,2));


%% get preference filter
minInfo = 1;
postConPref = postPrefMinusPre(:,:,1)./sum(postPrefMinusPre,3);
preConPref = prePrefCon(:,:,1)./sum(prePrefCon,3);
minPref = min(max(postPrefMinusPre,[],3),max(prePrefCon,[],3));

isPrefInfo = ~isnan(postConPref) & ~isnan(preConPref);
%isPrefInfo = isPrefInfo & ~isnan(synPref) & ~isinf(synPref);
%isPrefInfo = isPrefInfo & (minPref>=minInfo);% & (predSyn>=minPredictedSynapseNumber);
% 
% subplot(2,2,3)
% image(postConPref*100);
% subplot(2,2,4)
% image(preConPref*100);


%% Filter by favorite
%%Produce binary masks for con matrix that filters according to seed
%%preference

prefThresh = .5;
sameSeedOneMat = (postConPref < prefThresh) & (preConPref < prefThresh) & isPrefInfo;
sameSeedTwoMat = (postConPref > prefThresh) & (preConPref > prefThresh) & isPrefInfo;
difSeedOneMat = (postConPref > prefThresh) & (preConPref < prefThresh) & isPrefInfo;
difSeedTwoMat = (postConPref < prefThresh) & (preConPref > prefThresh) & isPrefInfo;

isSame = sameSeedOneMat | sameSeedTwoMat;
isDif = difSeedOneMat | difSeedTwoMat;
% 
% sameSeedCon = con(isSame);
% difSeedCon = con(isDif);

%% Collect masks

mask.isPrefInfo = isPrefInfo>0;

tracedPost = repmat(sum(predSyn,1),[size(sourceCon,1) 1]);

mask.tracePost = tracedPost>0;
mask.sameSeedOneMat = sameSeedOneMat>0;
mask.sameSeedTwoMat = sameSeedTwoMat>0;
mask.isSame = isSame>0;
mask.isDif = isDif>0;



%% Get predictions
binRes = 5;
%predSyn = overlapPredictSynOfSubset(skelOverlapPred,binRes,(mask.isPrefInfo & mask.tracePost));
%checkSynError = conError(con,predSyn,(mask.isPrefInfo & mask.tracePost))
synPref = (sourceCon )./predSyn;

mask.predOne = predSyn>=1;
% image(predSyn*20)
% image(mask.predOne*1000)

%% Group preferences
sameSeedPred = predSyn(mask.isPrefInfo & mask.tracePost & isSame);
difSeedPred = predSyn(mask.isPrefInfo & mask.tracePost & isDif);
sameSeedPref = synPref(mask.isPrefInfo & mask.tracePost & isSame & mask.predOne);
difSeedPref = synPref(mask.isPrefInfo & mask.tracePost & isDif & mask.predOne);

%% preference stats
medPref = [mean(sameSeedPref)
mean(difSeedPref)];

if 1
%% bar preferences
% clf

conBin = [0:1:10];
%conBin = [min([sameSeedPref; difSeedPref]): .5 :max([sameSeedPref; difSeedPref])+1];
histSame = histc(sameSeedPref,conBin)/length(sameSeedPref) * 100;
histDif = histc(difSeedPref,conBin)/length(difSeedPref)*100;
% 
% bar(conBin(1:end),[histSame(1:end) histDif(1:end)],'barWidth',1)
% 
% mean(sameSeedPref),std(sameSeedPref)/sqrt(length(sameSeedPref))
% mean(difSeedPref),std(difSeedPref)/sqrt(length(difSeedPref))
% ranksum(sameSeedPref,difSeedPref)
% 
% median(sameSeedPref)
% rangeX(sameSeedPref)
% median(difSeedPref)
% rangeX(difSeedPref)

end


prefStat.medPref = medPref;
prefStat.histSame = histSame;
prefStat.histDif = histDif;
prefStat.Ndif = length(difSeedPref);
prefStat.Nsame = length(sameSeedPref);
prefStat.N = prefStat.Ndif  + prefStat.Nsame ;
