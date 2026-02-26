function[diSub] = dilateSub(sub,ball)


%%
ballInd = find(ball);
[sY sX sZ] = size(ball);
maxBall = max([sY sX sZ]);
[bY bX bZ] = ind2sub(size(ball),ballInd);

%sub = subRaw+maxBall; %% buffer subs to stay above 0
subSize = max(sub,[],1)+maxBall;
allInds = zeros(prod(subSize),1,'uint8');

%%
% maxSize = 1000000;
%
% groupSize = floor(maxSize/length(bY));
groupSize = 1000; %group subs
shiftR = 1:(size(sub,1)/groupSize+1);

shiftInds = zeros(groupSize,length(bY));
shiftSub = zeros(groupSize,3,length(bY));
uniqueInds = [];

for r = 1:length(shiftR)
    %disp(sprintf('%d of %d',r,length(shiftR)))
    iRange = (shiftR(r)-1)*groupSize + [1:groupSize];
    iRange(iRange>size(sub,1)) = size(sub,1);
    subRange = sub(iRange,:);
    
    parfor i = 1:length(bY)
        shiftSub(:,:,i) = [subRange(:,1)+bY(i) subRange(:,2)+bX(i) subRange(:,3) + bZ(i)];
        
    
    end
    shiftInds = sub2ind(subSize,shiftSub(:,1,:),shiftSub(:,2,:),shiftSub(:,3,:));
    allInds(shiftInds(:)) = 1;
    
end

uniqueInds = find(allInds);

%%

vol1 = zeros(subSize,'uint8');
vol1(uniqueInds) = 1;
sumVol = squeeze(sum(vol1,1));
image(sumVol*50/size(ball,1))
%%
[dY dX dZ] = ind2sub(subSize,uniqueInds);
diSub = [dY dX dZ];
diSub = diSub - (maxBall+1)/2;

diSub(diSub<1) = 1;












