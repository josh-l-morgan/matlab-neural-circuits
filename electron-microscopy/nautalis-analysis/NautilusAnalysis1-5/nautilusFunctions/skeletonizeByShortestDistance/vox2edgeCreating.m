
%% Test Data

obj = ones(10,10,10);
inds = find(obj>0);
[y x z] = ind2sub(size(obj),inds);
subs = [y x z];

%% study
vNum = length(inds);


%% thourough con1
con1 = zeros(vNum,6);
conThresh = 1.9;

for i = 1:vNum;
    dists = sqrt(sum((subs-repmat(subs(i,:),[vNum 1])).^2,2));
    conTo = find((dists >0) & (dists<conThresh));
    con1(i,1:length(conTo)) = conTo;
    dist1(i,1:length(conTo)) = dists(conTo);
end


%% shortest path
con2 = con1;
seed = 1;
dist2seed = ones(1,vNum)*vNum;
dist2seed(seed) = 0;
pred = zeros(vNum,1);

nextCon = seed; %initialize next verticies to try

for vRep = 1:vNum
    lastCon = nextCon;
    nextCon = [];
    for i = 1:length(lastCon) % run all next verticies
        cV = lastCon(i); % grab current vertex
        conTo = con2(cV,:); % find who is connected to current vertex
        con2(cV,:) = 0; % clear connections once vertex is run once
        conTo = conTo(conTo>0); % grab actual connections
        distTo = dist1(cV,conTo>0); %
        betterPath = ( dist2seed(conTo) > (dist2seed(cV)+distTo)); %find connections to which better path has been found
        dist2seed(conTo(betterPath)) = (dist2seed(cV)+1); % Change distances to improve
        pred(conTo(betterPath)) = cV; % mark improved position by predecessor.
        
        nextCon = cat(2,nextCon,conTo); %  add connected to list for next start verticies.
    end
    
    if ~sum(con2(:,1))
        break
    end
end


showDisplay = 1;
if showDisplay

%% Draw paths56)
colormap gray(256)
for p = 1:vNum
showPath = obj*0;

tip = p;

nextPred = tip;
for v = 1:vNum
    if ~nextPred
        break
    end
    if ~showPath(nextPred)
        showPath(nextPred) = 1;
    end
    nextPred = pred(nextPred);
end
image(sum(showPath,3)*30);
pause(.1)

end
    
end % if display

