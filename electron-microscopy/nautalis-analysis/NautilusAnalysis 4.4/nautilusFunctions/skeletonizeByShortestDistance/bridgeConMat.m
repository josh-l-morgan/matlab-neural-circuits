function[allVox] = bridgeConMat(allVox)


clf
subs = dropSubs(allVox.subs);
showSteps = 1;

seed = allVox.seedInd;

renderConK(subs,[],[1 0 0])

allVox.conMatOld = allVox.conMat;
bSubs = dropSubs(allVox.subs);
rangeA = max(bSubs,[],1);
allInds = sub2ind(rangeA,bSubs(:,1),bSubs(:,2),bSubs(:,3));

con = allVox.conMat;

segID = zeros(size(con,1),1);
lastSegID = segID;

segCount = 0;
pairs = [];
clf
while 1
    segCount = segCount + 1;
    segID(seed) = segCount;
    
    %%Spread
    while 1
        for d = 1:26
            moveID =  (segID>0) & (con(:,d)>0);
            segID(con(moveID,d)) = segID(moveID);
        end
        %         image(segID*100)
        %         pause(.01)
        if sum(abs(lastSegID-segID))==0
            break
        end
        lastSegID =segID;
    end
    
    %%Find nearest
    isFilled = segID>0;
    isNot =segID==0;
    fullPosF = subs(isFilled,:);
    fullPosN = subs(isNot,:);
    clf
    if showSteps, renderConK(fullPosF,[],[1 1 1],.1), end
  
    
    if ~sum(isNot)
        break
    end
    if showSteps, renderConK(fullPosN,[],[0 0 1],.1),pause(.01), end
    
    
    
    %% Find potential overlaps
    
    boxSize = 200;
    rangeF = [min(fullPosF,[],1)-10;max(fullPosF,[],1)+10];
    rangeN = [min(fullPosN,[],1)-10; max(fullPosN,[],1)+10];
    rangeBuf = 1;
    
    for b = 1:20
        
        closeN = [];
        closeF = [];
        for y = 1:ceil(rangeN(2,1)/boxSize)
            for x = 1:ceil(rangeN(2,2)/boxSize)
                for z = 1:ceil(rangeN(2,3)/boxSize)
                    inRange =  ((subs(:,1)+rangeBuf)>((y-1)*boxSize)) & ((subs(:,1)-rangeBuf)<(y*boxSize)) & ...
                        ((subs(:,2)+rangeBuf)>((x-1)*boxSize)) & ((subs(:,2)-rangeBuf)<(x*boxSize)) & ...
                        ((subs(:,3)+rangeBuf)>((z-1)*boxSize)) & ((subs(:,3)-rangeBuf)<(z*boxSize));
                    
                    posF = subs(isFilled & inRange,:);
                    posN = subs(isNot & inRange,:);
                    
                    if ~isempty(posF) & ~isempty(posN)
                        
                        subRangeF = [min(posF,[],1)-10;max(posF,[],1)+10];
                        subRangeN = [min(posN,[],1)-10; max(posN,[],1)+10];
                        
                        
                        inRangeF =  ((posF(:,1)+rangeBuf)>((subRangeN(1,1)))) & ((posF(:,1)-rangeBuf)<(subRangeN(2,1))) & ...
                            ((posF(:,2)+rangeBuf)>((subRangeN(1,2)))) & ((posF(:,2)-rangeBuf)<(subRangeN(2,2))) & ...
                            ((posF(:,3)+rangeBuf)>((subRangeN(1,3)-1))) & ((posF(:,3)-rangeBuf)<(subRangeN(2,3)));
                        
                        inRangeN =  ((posN(:,1)+rangeBuf)>((subRangeF(1,1)))) & ((posN(:,1)-rangeBuf)<(subRangeF(2,1))) & ...
                            ((posN(:,2)+rangeBuf)>((subRangeF(1,2)))) & ((posN(:,2)-rangeBuf)<(subRangeF(2,2))) & ...
                            ((posN(:,3)+rangeBuf)>((subRangeF(1,3)-1))) & ((posN(:,3)-rangeBuf)<(subRangeF(2,3)));
                        
                        closeF = cat(1,closeF,posF(inRangeF,:));
                        closeN = cat(1,closeN,posN(inRangeN,:));
                        
                    end %is pos
                    
                end %run z
            end %run x
        end %runy
        
        if ~isempty(closeF) & ~isempty(closeN)
            
            uF = unique(sub2ind(rangeA,closeF(:,1),closeF(:,2),closeF(:,3)));
            uN = unique(sub2ind(rangeA,closeN(:,1),closeN(:,2),closeN(:,3)));
            [y1 x1 z1] = ind2sub(rangeA,uF);
            closeF = [y1 x1 z1];
            [y2 x2 z2] = ind2sub(rangeA,uN);
            closeN = [y2 x2 z2];
            if showSteps, renderCon(closeF+.5,[],[1 0 0]), end
            if showSteps, renderCon(closeN+.5,[],[0 1 0]), end
            pause(.01)
            break
            
        end
        
        rangeBuf = rangeBuf * 2;
        
    end
    
    
    %% 
       
    %dists = sqrt((posF(:,1)-posN(:,1)').^2 + (posF(:,2)-posN(:,2)').^2 + (posF(:,3)-posN(:,3)').^2);
    dists = sqrt((closeF(:,1)-closeN(:,1)').^2 + (closeF(:,2)-closeN(:,2)').^2 + (closeF(:,3)-closeN(:,3)').^2);
    
    
    
    minDist = min(dists(:));
    [f n] = find(dists == minDist,1);
    
    hitF = closeF(f,:);
    hitN = closeN(n,:);
    hitFind = sub2ind(rangeA,hitF(1),hitF(2),hitF(3));
    hitNind = sub2ind(rangeA,hitN(1),hitN(2),hitN(3));
       
    pF = find((allInds==hitFind) & isFilled,1);
    renderCon(subs(pF,:)+.5,[],[1 1 0],1)
    pN =  find((allInds==hitNind) & ~isFilled);
    renderCon(uniqueSubs(subs(pN,:))+.5,[],[0 1 1],1)
    pause(.01);
    pN = setdiff(pN,pF); % in case of multiple hits
    pN = pN(1);
    pairs  = cat(1,pairs,[pF pN]);
    
    seed = pN; %ready for next spreading
    
    
    if showSteps,
        hold on
        for b = 1:size(pairs,1)
            nList = pairs(b,:);
            runE = bSubs(nList,:);
            plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','w')
        end
        
        pause(.01)
    end
    
end

%% change conMat

for p = 1:size(pairs,1)
    
    con1 = con(pairs(p,1),:);
    con2 = con(pairs(p,2),:);
    
    empty1 = find(con1==0,1);
    empty2 = find(con2==0,1);
    
    con(pairs(p,1),empty1) = pairs(p,2);
    con(pairs(p,2),empty2) = pairs(p,1);
    
end

allVox.conMat = con;
allVox.voxBridges = pairs; 

if 1 %% show
    clf
    renderConK(subs,allVox.conMat,[1 0 0])
    hold on
        for b = 1:size(pairs,1)
            nList = pairs(b,:);
            runE = bSubs(nList,:);
            plot3(runE(:,1),runE(:,2),runE(:,3),'linewidth',3,'color','w')
        end  
    
    ax = gca;
    ax.Clipping = 'off';
    
end

pause(.01)








