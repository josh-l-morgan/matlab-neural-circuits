

%generic connectivity skew test
%Josh Morgan

%% options
reps = 10000;

%% make dummy data until someone gives me real data
preIn = fix(rand(100,1)*90) + 1;
postIn = fix(rand(100,1)*3)+1;

%% renumber IDs
ids = sort(unique(preIn));
newIds = 1:length(ids);
lookupIds(ids) = newIds;
pre = lookupIds(preIn);

ids = sort(unique(postIn));
newIds = 1:length(ids);
lookupIds(ids) = newIds;
post = lookupIds(postIn);


%% get data info
preNum = length(unique(pre));
postNum = length(unique(post));
synNum = length(pre);

postDist = hist(post,postNum)/synNum;
preDist = hist(pre,preNum)/synNum;

%% get ax positions
for a = 1:preNum;
    axPos{a} = find(pre == a);
end

%% run 
for i = 1:reps
    pos = rand(synNum,1);
    [a newPos] = sort(pos);
    newDend = post(newPos);
    for a = 1:preNum
        hitPost = newDend(axPos{a});
        if length(hitPost)>1
            aHist = hist(hitPost,[.5:1:(postNum-.5)]);
            difD = (max(aHist)-min(aHist) -1);
            difD(difD<0) = 0;
            
        else 
           difD = 0; % no skew
        end
       
        difs(a) = difD;
    end
    mean(difs);
    aveSkew(i) = mean(difs);
    
    if mod(i,1000) == 0
        sprintf('running %d of %d reps',i,reps)
    end
end

skewBins = [min(aveSkew):(max(aveSkew)-min(aveSkew))/100:max(aveSkew)];
hist(aveSkew,skewBins);
bar(rDist)
    

