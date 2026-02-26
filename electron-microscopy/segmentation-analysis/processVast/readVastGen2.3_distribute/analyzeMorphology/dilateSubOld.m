function[diSub] = dilateSub(subRaw,ball)


%%
ballInd = find(ball);
[sY sX sZ] = size(ball);
maxBall = max([sY sX sZ]);
[bY bX bZ] = ind2sub(size(ball),ballInd);

sub = subRaw+maxBall;
subSize = max(sub,[],1)+maxBall;

%%
% maxSize = 1000000;
% 
% groupSize = floor(maxSize/length(bY));
groupSize = 1000;
shiftR = 1:(length(bY)/groupSize+1);


shiftInds = zeros(size(sub,1),groupSize);
uniqueInds = [];
for r = 1:length(shiftR)
    %disp(sprintf('%d of %d',r,length(shiftR)))
    iRange = (shiftR(r)-1)*100 + [1:groupSize];
    iRange = iRange(iRange<=length(bY));
    sbY = bY(iRange); sbX = bX(iRange); sbZ = bZ(iRange);
    
parfor i = 1:length(iRange)
   shiftSub = [sub(:,1)+sbY(i) sub(:,2)+sbX(i) sub(:,3) + sbZ(i)];
   shiftInds(:,i) = sub2ind(subSize,shiftSub(:,1),shiftSub(:,2),shiftSub(:,3));
   %allInds =  unique([allInds; shiftInds]);
end

uniqueInds = unique(cat(1,uniqueInds,shiftInds(:)));


%uniqueInds = cat(1,uniqueInds, unique(shiftInds(:)));
%uniqueInds = [uniqueInds; shiftInds(:)];

end

uniqueInds = uniqueInds(uniqueInds>0);

%%

vol1 = zeros(subSize,'uint8');
vol1(uniqueInds) = 1;
sumVol = squeeze(sum(vol1,1));
image(sumVol*50/size(ball,1))
%%
[dY dX dZ] = ind2sub(subSize,uniqueInds);
diSub = [dY dX dZ];
diSub = diSub - maxBall;

diSub(diSub<1) = 1;












