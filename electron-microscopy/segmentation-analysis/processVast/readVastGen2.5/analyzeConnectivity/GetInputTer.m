
clear all
MPN = GetMyDir
TPN = [MPN 'manyTerSubs\'];

downSamps = [ 1  1 4 8  64]
lookDists = [ 0  4 4  8  16]

downSamps = [  1 ]
lookDists = [  0 ]

lookMicrons =  .03 * 4  + lookDists .* downSamps * 8 * 0.03

subplot(1,1,1)
for s = 1: length(lookDists)
    
downSamp = downSamps(s);
lookDist = lookDists(s);
fileName = sprintf('terSubs_Ds%d_Ds%d_Look%d.mat',8,downSamp,lookDist)

load([TPN fileName]);
synMat = terSubs.synMat;
touchMat = terSubs.touchMat;

preSubs = [];
for i = 1: length(terSubs.preList)
    subs = terSubs.cells(terSubs.cellList == terSubs.preList(i)).terSubs;
	preSubs = cat(1,preSubs,subs);
end

postSubs = [];
for i = 1 : length(terSubs.postList)
   subs = terSubs.cells(terSubs.cellList == terSubs.postList(i)).terSubs;
    postSubs = cat(1,postSubs,subs);
end

maxSubs = max(cat(1,postSubs,preSubs),[],1);
preInds = sub2ind(maxSubs,preSubs(:,1),preSubs(:,2),preSubs(:,3));
postInds = sub2ind(maxSubs,postSubs(:,1),postSubs(:,2),postSubs(:,3));

allPrePostSubs{1} = preSubs;
allPrePostSubs{2} = postSubs;

%% Get overlap inds
clear sameSubs
for i = 1 : length(terSubs.postList)
   subs = terSubs.cells(terSubs.cellList == terSubs.postList(i)).terSubs;
   postInd = sub2ind(maxSubs,subs(:,1),subs(:,2),subs(:,3));
   sameInds{i} = intersect(preInds,postInd);
   [y x z] = ind2sub(maxSubs,sameInds{i});
   sameSubs{i} = [y x z];
end


%%
for i = 1:length(sameSubs)
   Isum = sumSub(sameSubs{i},1,maxSubs); 
   Isum = imresize(Isum,.1)*100;
   image(Isum*100)
   pause(.1)
    
end

    Ifull = showSubSums(sameSubs);
    image(Ifull*1.4)

%%
obNum = length(sameSubs)
colMap = hsv(256);
col = colMap(ceil((1:obNum)*256/obNum),:);
col1 = col(randperm(size(col,1)),:);


obNum = length(allPrePostSubs)
colMap = hsv(256);
col = colMap(ceil((1:obNum)*256/obNum),:);
col2 = col(randperm(size(col,1)),:);



allSubs = cat(1,sameSubs{:});
viewR = 30;
steps = min(allSubs(:,1)):5:max(allSubs(:,1));


for i =1:length(steps)
    mid = steps(i);
    sr = [ mid-viewR 0 0; mid+viewR maxSubs(2) maxSubs(3) ];
    Ifull = showSubSums(sameSubs,sr,col1);
    Ipost = showSubSums(allPrePostSubs,sr,col2);
    Iboth = Ifull*1.7 + Ipost*.3;
    Iboth = imresize(Iboth,.5,'bicubic');
    image(Iboth)
    if i == 1;
        recBoth = zeros(size(Iboth,1),size(Iboth,2),3,length(steps),'uint8');
    end
    recBoth(:,:,:,i) = Iboth;
    pause(.01)
end
%%

while 1
for i = 1:size(recBoth,4)
   image(recBoth(:,:,:,i)),
   pause(.08) 
end
for i = size(recBoth,4):-1:1
   image(recBoth(:,:,:,i)),
   pause(.08) 
end
end



%% Get overlap inds
clear sameSubs
for i = 1 : length(terSubs.preList)
   subs = terSubs.cells(terSubs.cellList == terSubs.preList(i)).terSubs;
   preInd = sub2ind(maxSubs,subs(:,1),subs(:,2),subs(:,3));
   sameInds{i} = intersect(preInd,postInds);
   [y x z] = ind2sub(maxSubs,sameInds{i});
   sameSubs{i} = [y x z];
end

%%

Ifull = showSubSums(sameSubs,sr);


obNum = length(sameSubs)
colMap = hsv(256);
col = colMap(ceil((1:obNum)*256/obNum),:);
col1 = col(randperm(size(col,1)),:);


obNum = length(allPrePostSubs)
colMap = hsv(256);
col = colMap(ceil((1:obNum)*256/obNum),:);
col2 = col(randperm(size(col,1)),:);

allSubs = cat(1,sameSubs{:});
viewR = 30;
steps = min(allSubs(:,1)):5:max(allSubs(:,1));


for i =1:length(steps)
    mid = steps(i);
    sr = [ mid-viewR 0 0; mid+viewR maxSubs(2) maxSubs(3) ];
    Ifull = showSubSums(sameSubs,sr,col1);
    Ipost = showSubSums(allPrePostSubs,sr,col2);
    Iboth = Ifull*1.7 + Ipost*.3;
    Iboth = imresize(Iboth,.5,'bicubic');
    image(Iboth)
    if i == 1;
        recBoth = zeros(size(Iboth,1),size(Iboth,2),3,length(steps),'uint8');
    end
    recBoth(:,:,:,i) = Iboth;
    pause(.01)
end
%%

while 1
for i = 1:size(recBoth,4)
   image(recBoth(:,:,:,i)),
   pause(.08) 
end
for i = size(recBoth,4):-1:1
   image(recBoth(:,:,:,i)),
   pause(.08) 
end
end








end